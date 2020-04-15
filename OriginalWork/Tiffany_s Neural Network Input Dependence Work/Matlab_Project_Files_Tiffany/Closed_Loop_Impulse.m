%% Tiffany Hamstreet - Portland State University
%% Animatlab Serial interface with MATLAB to run Inverted Pendulum Simulation (credit to Wade original HEBI/Matlab/AM code)
%% For IMPULSE RESPONSE with CLOSED LOOP PD Controller/Inverted Pendulum Simulation, includes serial connections to Error (N3), Kp Scope, and Kd Scope Neurons

clear all;
clc;

% setup serial connection properties
delete(instrfindall)
r = serial('COM2');     % receive data from AM over COM1
set(r,'BaudRate',256000); % it is critical that strict baudrate emulation is enabled on the virtual serial port driver
r.Timeout = 50; % default 10
fopen(r);
s = serial('COM3');     % send data to AM over COM2
set(s,'BaudRate',256000); % Animatlab seems to be tolerant of very high virtual baudrates. I could have gone to 10^6
fopen(s);

% initialize variables for serial communication
n = 200000; % initialize the maximum number of serial bytes that will be read before vectors are full (set to well above what sim time needs)

% initialize simulation variables and parameters
Sample_Frequency = 200;  % sample frequency Hz
dt = 1/Sample_Frequency; %s
% rad_to_nA = 1.5*(6/pi()); % map radians to current, 30 deg = pi/6 = 1.5 nA
% nA_to_rad = pi()/(1.5*6);
convert_to_nA = 1000;
loopTime = zeros(n,1);
% gain = 833;  % amplifier gain out of controller to simulation (weird offset happens at gain = 200), higher gain, faster settling time
gain = 1388;
Veq = -.060; %V used to normalize commanded torque

pos_head = [255 255 1 18 0 23 0]; %defines header, message ID (always 1), message size (2 bytes) and Data ID for Actual Angle = 23 in AM (2 bytes), starting actual angle = 0
pos_des_head = [24 0]; %[data ID for desired angle = 24 in Animatlab Input Serial IO properties, desired position = 0]
pos_msg = zeros(1,12);
crit_error = 0;

% auto bode frequencies and input amplitudes to loop through
theta_impulse = [1];  % IMPULSE amplitude in nA
N_impulses = 1;

delay = 0.1;              % s
delaySteps = delay/dt;      % delay must be divisible into integers by dt
% % set up loop timer (robotics toolbox)
% desiredRate = Ts;
% rateObj = robotics.Rate(desiredRate);
% rateObj.OverrunAction  = 'slip';

% reset(rateObj);

% Auto Bode settings
n_columns = 1;
    
%% Pull Torques from AM, run IP simulation, Send Actual Position to AM

for bb=1:N_impulses
    
%     NNoutput = zeros(N_test_cycles*N_points_per_cycle,N_frequencies);    % sum of NN output
%     NNinput_des = zeros(N_test_cycles*N_points_per_cycle,N_frequencies); % 
%     time = zeros(N_test_cycles*N_points_per_cycle+1,N_frequencies);

    NNoutput = zeros(n,n_columns);      % sum of NN output
    NNinput_des = zeros(n,n_columns);   %
    NNinput_angle = zeros(n,n_columns); % position as current to apply to AM (zero to N1)
%     NNinput_angle(50,1) = 1;
%     NNinput_des(50,1) = 1;             % position as 1nA current to apply impulse
%     NNinput_des(921,1) = -2;
    time = zeros(n,n_columns);
    CW_T = zeros(n,n_columns);          % CW torque commmand vector
    CCW_T = zeros(n,n_columns);         % CCW torque command vector
    Error_N3 = zeros(n,n_columns);      % Error Neuron 3 membrane voltage
    KpError_N10 = zeros(n,n_columns);   % Kp*Error Neuron 10 membrane voltage
    KddError_N13 = zeros(n,n_columns);  % Kp*dError Neuron 13 membrane voltage
    theta = zeros(n,2); % position (radians?) and velocity (rad/s?)
    theta(2,1) = theta_impulse;                   % initial displacement ~2.8 degrees
    Tapp = zeros(n,1);  % Torque applied (controller output*gain)
    
    % Default starting values at first index
    time(2,:) = dt;
    CW_T(1,:) = Veq;
    CCW_T(1,:) = Veq;
    Error_N3(1,:) = Veq;
    KpError_N10(1,:) = Veq;
    KddError_N13(1,:) = Veq;
    
    tend = 50;      % simulation time
    rowIndex = 2;

%         fprintf('start frequency %i Hz \n',frequencies_Hz(jj));       
    jj = 1;

    while time(rowIndex,jj) < tend
        tic
            
            Data_Vec = [CCW_T(rowIndex-1,jj) ...
                CW_T(rowIndex-1,jj) ... 
                Error_N3(rowIndex-1,jj) ... 
                KddError_N13(rowIndex-1,jj) ...
                KpError_N10(rowIndex-1,jj)];

            ID_Vec = [21, 22, 25, 26, 27];  % [CCW_T CW_T Error KddError KpError]

            New_Data = GetAnimatData(r,Data_Vec,ID_Vec,rowIndex);

            CCW_T(rowIndex,jj) = New_Data(1);
            CW_T(rowIndex,jj) = New_Data(2);
            Error_N3(rowIndex,jj) = New_Data(3);
            KddError_N13(rowIndex,jj) = New_Data(4);
            KpError_N10(rowIndex,jj) = New_Data(5);   

            NNoutput(rowIndex,jj) = ((CCW_T(rowIndex,jj)-Veq)-(CW_T(rowIndex,jj)-Veq))*convert_to_nA; % total control torque from AM, adjusted by equilibrium voltage multiplied by amplifier gain

                if abs(NNoutput(rowIndex,jj) - NNoutput(rowIndex-1,jj)) > 5
                    fprintf('Serial Error in Frequency %d, rowIndex %d \n',jj,rowIndex)
                end

                if rowIndex > delaySteps
                    Tapp(rowIndex) = NNoutput(rowIndex-delaySteps,jj)*gain; %convert theta_out (input to controller to current
                end

                % Inverted Pendulum Simulation
                if rowIndex > 2
                    [theta(rowIndex,:)] = Simulate_InvertedPend_Tiff(theta(rowIndex-1,:),Tapp(rowIndex),dt,dt,1); %[position, velocity] inputs(current position/velocity, current corrective Torque applied, simulation time, dt, howmuch) for howmuch = 1, simulation returns only final calculated value, so with time = dt, simulation only runs through one time each time it's accessed in the main code
                end
%               
                NNinput_angle(rowIndex,jj) = theta(rowIndex,1); %convert theta_out (input to controller to current
%               
            % write new desired position to Animatlab
            pos_byte = typecast(single(NNinput_angle(rowIndex,jj)),'uint8');  %loads position current to send to AM
%                 des_byte = typecast(single(pos_des(cmd+1)),'uint8'); %loads desired position to send to AM
            des_byte = typecast(single(NNinput_des(rowIndex,jj)),'uint8'); %loads desired position to send to AM
            chksum = mod(sum([pos_head pos_byte pos_des_head des_byte]),256);
            pos_msg = [pos_head pos_byte pos_des_head des_byte chksum]; % package up serial message
            fwrite(s, pos_msg, 'uint8'); % write serial message to animatlab sim 

        %flush serial buffers
        flushinput(r);
        flushoutput(s);
        time(rowIndex+1,jj) = time(rowIndex,jj) + dt;   % ms           
        rowIndex = rowIndex + 1;

        % if not using Robotics toolbox, can use this tic toc timing to delay
        % loop until end of timestep
        z = toc;    

        if z < dt
            pause(dt - z);
            z = 0;
        else
            z = 0;
        end

    %     waitfor(rateObj);
        loopTime(rowIndex,1) = toc;

    end

    Error_N3(:,jj) = (Error_N3(:,jj) - Veq)*convert_to_nA;
    KpError_N10(:,jj) = (KpError_N10(:,jj) - Veq)*convert_to_nA;
    KddError_N13(:,jj) = (KddError_N13(:,jj) - Veq)*convert_to_nA;

    timeNNinput =       [time(1:50000,1)    NNinput_des(1:50000,:)  ];  % desired position input in nA
    timeError_N3 =      [time(1:50000,1)    Error_N3(1:50000,:)     ];  % normalized to equ voltage in nA
    timeKpError_N10 =   [time(1:50000,1)    KpError_N10(1:50000,:)  ];  % normalized to equ voltage in nA
    timeKddError_N13 =  [time(1:50000,1)    KddError_N13(1:50000,:) ];  % normalized to equ voltage in nA
    timeNNoutput =      [time(1:50000,1)    NNoutput(1:50000,:)     ];  % output (proportional to torque but pre-gain) in nA
    timeNNangleOut =    [time(1:50000,1)    NNinput_angle(1:50000,:)];  % actual position in nA

    % UPDATE FORMAT SPEC TO #FREQ + 1 AND UPDATE FILENAME
    formatSpec = '%f %f\n'; % for 1 input frequency
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ImpCL_Outputs\NNinput_Test(',num2str(theta_impulse(bb)),').dat'],'w');
    [row,c] = size(timeNNinput);
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeNNinput(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ImpCL_Outputs\ErrorOutput_Test(',num2str(theta_impulse(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeError_N3(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ImpCL_Outputs\KpOutput_Test(',num2str(theta_impulse(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeKpError_N10(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ImpCL_Outputs\KdOutput_Test(',num2str(theta_impulse(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeKddError_N13(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ImpCL_Outputs\NNOutput_Test(',num2str(theta_impulse(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeNNoutput(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ImpCL_Outputs\NNangleOut_Test(',num2str(theta_impulse(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeNNangleOut(zz,:));
    end
    
    fclose(fid(bb));
    
end

fprintf('simulation done \n');

% loopData = statistics(rateObj)
% delete(rateObj)

%% time domain plots for sanity checks

% figure
% plot(time(:,1),theta(:,1))
% title('IP Angle Theta in rads')
% xlim([0 15])
% ylim([-1 1])
% xlabel('time')
% ylabel('radians')

figure
plot(time(:,1),Error_N3(:,1))
title('Error Signal')
xlim([0 60])
ylim([-1 1])
xlabel('time')
ylabel('nA')

figure
plot(time(:,1),KpError_N10(:,1))
title('Kp*Error Signal')
xlim([0 60])
ylim([-1 1])
xlabel('time')
ylabel('nA')

figure
plot(time(:,1),KddError_N13(:,1))
title('Kd*dError Signal')
xlim([0 60])
ylim([-1 1])
xlabel('time')
ylabel('nA')

figure
plot(time(:,1),NNoutput(:,1))
title('NN Output Signal')
xlim([0 60])
ylim([-1 1])
xlabel('time')
ylabel('nA')

figure
hold on
plot(time(:,1),NNinput_des(:,1))
plot(time(:,1),NNinput_angle(:,1))
legend('Command Angle (in nA)','Actual Angle (in nA)')
xlim([0 60])
ylim([-1 1])
xlabel('time')
ylabel('nA')

%% FFT(impulse)

N = length(NNinput_angle(:,1));

% Peterka's "normal subject" constants
m = 83.3; %kg 
h = .896; %m - center of mass height
I = 81.1; %kg m^2
g = 9.81;% m/s

% closed loop
m1 = floor(N/2)+1;  % length of post-DFT signal
z_fft = fft(NNinput_angle(:,1));
p = abs(z_fft);  %two sided spectrum
fRes_pre = p(1:m1);   %single-sided spectrum
fRes = fRes_pre/N;    %scale by 1/N
fRes(2:end-1) = 2*fRes(2:end-1);    %scale all but first and last value by 2
phase = angle(z_fft);
f = linspace(0,Sample_Frequency/2,m1);

s = tf('s');

Kp = 0.21;
Kd = 0.050;
tau = 0.15;
gain = 5000;                % amplifier gain applied to control signal
[num,den] = pade(tau,1);    % pade delay approximation first order (in order to use impulse command)
delay = tf(num,den);
Gc = Kp + Kd*s;           % controller TF
Gp = 1/(I*s^2-m*g*h);       % Inverted Pendulum TF
oltf = delay*Gc*Gp*gain;
cltf = feedback(oltf,1);

% bode data for comparison
[magCLTF,phaseCLTF,freqCLTF] = bode(cltf);
magCLTF = squeeze(magCLTF);
phaseCLTF = squeeze(phaseCLTF);
freqCLTF = freqCLTF/(2*pi);

% open loop
z2_fft = fft(NNoutput(:,1));
p = abs(z2_fft);  %two sided spectrum
fRes2_pre = p(1:m1);   %single-sided spectrum
fRes2 = fRes2_pre/N;    %scale by 1/N
fRes2(2:end-1) = 2*fRes2(2:end-1);    %scale all but first and last value by 2
phase2 = angle(z2_fft);

% Closed loop Frequency Response
figure
semilogx(f,fRes_pre)  % used pre-scaled spectral density
hold on
semilogx(freqCLTF,magCLTF)
% plot(freqCLTF,magCLTF,'Linewidth',1.2)
xlabel('Frequency (Hz)')
ylabel('Spectral Density')
title('Closed Loop PD/IP Frequency Response fft(Impulse Response) of Actual Position')
xlim([.01 10])

% %Open loop Frequency Response
% figure
% % hold on
% semilogx(f,fRes2_pre)  % used pre-scaled spectral density
% xlabel('Frequency (Hz)')
% % ylabel('Spectral Density')
% title('Open Loop PD Frequency Response fft(Impulse Response) of NN Output')
% xlim([.01 100])

%% STILL OPEN LOOP PLOTS HERE
% s = tf('s');
% 
% Kp = 1;
% Kd = 0.04;
% unknownGain = 1;
% tau = 0.004;
% 
% fb1 = 5;         % cutoff frequency in Hz (base 2)
% lowFiltDerivative = (1/(s/(2*pi*fb1)+1));  % low pass filter
% 
% fb2 = 80;         % cutoff frequency in Hz (base 25)
% lowFiltNeuron = (1/(s/(2*pi*fb2)+1));  % low pass filter
% 
% [num,den] = pade(tau,1);    % pade delay approximation first order (in order to use impulse command)
% delay = tf(num,den);
% 
% % sys = unknownGain*(Kp*lowFiltNeuron + Kd*s*lowFiltNeuron*lowFiltDerivative);
% sys2 = lowFiltNeuron^4*(Kp + Kd*s*lowFiltDerivative);
% sysDelay = delay*sys2;
% 
% [magDelay,phaseDelay] = bode(sysDelay,frequencies);    % requires frequency in rad/s
% magDelay = squeeze(magDelay);
% phaseDelay = squeeze(phaseDelay);
% phaseDelay = phaseDelay - 360;
% 
% phase_save_p(1) = phase_save_p(1) + 360;
% 
% % [mag,phase] = bode(tf1,frequencies);
% % mag = 20*log10(squeeze(mag));
% % phase = squeeze(phase);
% % 
% % opts = bodeoptions('cstprefs');
% % opts.FreqUnits = 'Hz';
% % figure
% % bode(sys*lowFiltDerivative*lowFiltNeuron,opts)
% % title('Baseline for Neural Controller: Open Loop PD Bode Kp = 1, Kd = 0.0319')
% % xlim([1 100])
% 
% figure
% semilogx(frequencies_Hz,Amplitude_save_p,'-o',frequencies_Hz,magDelay,'-x')
% title('Amplitude')
% xlabel('Frequencies Hz')
% ylabel('Mag Ratio (abs)')
% legend('Actual','Fit','Location','northwest')
% 
% figure 
% semilogx(frequencies_Hz,phase_save_p,'-o',frequencies_Hz,phaseDelay)
% title('Phase')
% xlabel('Frequencies Hz')
% ylabel('Phase degrees')
% % ylim([-100 10])
% 
% figure
% semilogx(frequencies_Hz,Amplitude_save_ratio_Kp,'m-+')
% title('Kp Out/Error In')
% xlabel('Frequencies Hz')
% ylabel('Mag Ratio (abs)')
% 
% figure
% semilogx(frequencies_Hz,Amplitude_save_ratio_Kd,'k-d')
% title('Kd*dE Out/Error In')
% xlabel('Frequencies Hz')
% ylabel('Mag Ratio (dB)')

%%
% figure
% semilogx(frequencies_Hz,Amplitude_in_save_Kp,frequencies_Hz,Amplitude_Out_save_Kp,'m-+')
% title('Kp Input & Output')
% xlabel('Frequencies Hz')
% ylabel('Amplitude (nA)')
% legend('Error Input Amplitude','Kp*Error Output Amplitude')
% 
% 
% figure
% semilogx(frequencies_Hz,Amplitude_in_save_Kd,frequencies_Hz,Amplitude_Out_save_Kd,'k-d')
% title('Kd Input & Output')
% xlabel('Frequencies Hz')
% ylabel('Amplitude (nA)')
% legend('Error Input Amplitude','Kd*dError Output Amplitude')

%% visualize the data

% Timer Statistics
% figure
% plot(loopTime(1:500),'o')
% title('50 Hz tic toc Run 5')
% ylabel('Sample time (s)')
% 
% mu = mean(loopTime(1:500))
% stdev = std(loopTime(1:500))
% coeffVar = stdev/mu