%% Tiffany Hamstreet - Portland State University
%% Animatlab Serial interface with MATLAB to run Inverted Pendulum Simulation (credit to Wade original HEBI/Matlab/AM code)
%% For automatic bode creation with CLOSED LOOP PD Controller/Inverted Pendulum Simulation, includes serial connections to Error (N3), Kp Scope, and Kd Scope Neurons

clear all;
clc;

% setup serial connection properties
delete(instrfindall)
r = serial('COM2');     % receive data from AM over COM1
set(r,'BaudRate',256000); % it is critical that strict baudrate emulation is enabled on the virtual serial port driver
r.Timeout = 45; % default 10
fopen(r);
s = serial('COM3');     % send data to AM over COM2
set(s,'BaudRate',256000); % Animatlab seems to be tolerant of very high virtual baudrates. I could have gone to 10^6
fopen(s);

% initialize variables for serial communication
n = 2000000; % initialize the maximum number of serial bytes that will be read before vectors are full (set to well above what sim time needs)

% initialize simulation variables and parameters
Veq = -.060; %V used to normalize commanded torque

Sample_Frequency = 200;     % sample frequency Hz
dt = 1/Sample_Frequency;    % s

delay = .10;                % s
delaySteps = delay/dt;      % delay must be divisible into integers by dt
% gain = 833;  % amplifier gain out of controller to simulation
gain = 1388;
% gain = 1500;
tend = 400;

simtime = 0:dt:tend-dt;
ntime = length(simtime)+1;
loopTime = zeros(ntime,1);
% rad_to_nA = 1.5*(6/pi()); % map radians to current, 30 deg = pi/6 = 1.5 nA

% setup serial package for writing to animatlab
pos_head = [255 255 1 18 0 23 0]; %defines header, message ID (always 1), message size (2 bytes) and Data ID for Actual Angle = 23 in AM (2 bytes), starting actual angle = 0
pos_des_head = [24 0]; %[data ID for desired angle = 24 in Animatlab Input Serial IO properties, desired position = 0]
pos_msg = zeros(1,12);

% Generate filtered Gaussian white noise for sys ID - length of signal
ranSig = 5*randn(ntime,1);  % the multiplier here is a tuned value to get mean ~ 1.4amp input (depends on sample and cutoff frequencies)

    % REDUCE CUTOFF FREQUENCY WHEN RESPONSE IS IN APPROPRIATE RANGE
    % Use lowpass Butterworth filter
    fc = 10;     % cutoff freq Hz
    order = 6;      % Butter filter order
    Wn = (fc/Sample_Frequency)/2; % cutoff frequency as ratio of sample freq

    [z,p,k] = butter(order,Wn,'low');
    [b,a] = zp2tf(z,p,k);
    [sos,g] = tf2sos(b,a);
%     freqz(sos,1024,Sample_Frequency)
    
    filtered_noise = filtfilt(sos,g,ranSig);
    
    time = 0:dt:tend;
%     x = zeros(length(time),1);
% 
%     f = linspace(0.01,2.5,100);
% 
%     for jj = 1:ntime
%         for ii = 1:length(f)
%             x(jj) = x(jj) + 0.1*sin(f(ii)*2*pi*time(jj));
%         end
%     end

NNinput_des = filtered_noise;
inputMean = mean(abs(NNinput_des))    

% figure
% plot(time,NNinput_des);

% auto bode frequencies and input amplitudes to loop through
% amp = [1 5 10];
amp = [round(mean(abs(NNinput_des)),1)];    % nA
N_amps = length(amp);
convert_to_nA = 1000;

N_frequencies = 1;

% % set up loop timer (robotics toolbox)
% desiredRate = Ts;
% rateObj = robotics.Rate(desiredRate);
% rateObj.OverrunAction  = 'slip';

% reset(rateObj);

%% Pull Torques from AM, run IP simulation, Send Actual Position to AM

for bb=1:N_amps

    NNoutput = zeros(ntime,N_frequencies);    % sum of NN output
%     NNinput_des = zeros(n,N_frequencies); %
    NNinput_angle = zeros(ntime,N_frequencies); % position as current to apply to AM (zero to N1)
    time = zeros(ntime,N_frequencies);
    CW_T = zeros(ntime,N_frequencies);          % CW torque commmand vector
    CCW_T = zeros(ntime,N_frequencies);         % CCW torque command vector
    Error_Scope = zeros(ntime,N_frequencies);      % Error Neuron 3 membrane voltage
    KpError_N10 = zeros(ntime,N_frequencies);   % Kp*Error Neuron 10 membrane voltage
    KddError_N13 = zeros(ntime,N_frequencies);  % Kp*dError Neuron 13 membrane voltage
    Error_N3 = zeros(ntime,N_frequencies);      % Error Neuron 3 membrane voltage
    Kp_scaler = zeros(ntime,N_frequencies);      % Kp scaler N8 membrane voltage
    Kd_scaler = zeros(ntime,N_frequencies);      % Kd scaler N11 membrane voltage
    Kp_syn_gain = zeros(ntime,N_frequencies);      % Kp syn gain membrane voltage
    Kd_syn_gain = zeros(ntime,N_frequencies);      % Kd syn gain membrane voltage
    kd_pre_gain = zeros(ntime,N_frequencies);      % kd*de/dt CW membrane voltage
    
    % Default starting values at first index
    time(2,:) = dt;
    CW_T(1,:) = Veq;
    CCW_T(1,:) = Veq;
    Error_Scope(1,:) = Veq;
    KpError_N10(1,:) = Veq;
    KddError_N13(1,:) = Veq;
    Error_N3(1,:) = Veq;
    Kp_scaler(1,:) = Veq;
    Kd_scaler(1,:) = Veq;
    kd_pre_gain(1,:) = Veq;
    
    for jj=1:N_frequencies

        theta = zeros(ntime,2); % position (radians) and velocity (rad/s)
        Tapp = zeros(ntime,1);  % Torque applied (controller output*gain)

        rowIndex = 2;
        
        fprintf('start simulation \n');       

        while time(rowIndex,jj) < tend
            tic

            Data_Vec = [CCW_T(rowIndex-1,jj) ...
                CW_T(rowIndex-1,jj) ... 
                Error_Scope(rowIndex-1,jj) ... 
                KddError_N13(rowIndex-1,jj) ...
                KpError_N10(rowIndex-1,jj) ...
                Error_N3(rowIndex-1,jj) ...
                Kp_scaler(rowIndex-1,jj) ...
                Kd_scaler(rowIndex-1,jj) ...
                Kp_syn_gain(rowIndex-1,jj) ...
                Kd_syn_gain(rowIndex-1,jj) ...
                kd_pre_gain(rowIndex-1,jj)];

            ID_Vec = [21, 22, 25, 26, 27, 28, 29, 30, 31, 32, 33];  % [CCW_T CW_T Error KddError KpError ...]

            New_Data = GetAnimatData(r,Data_Vec,ID_Vec,rowIndex);

            CCW_T(rowIndex,jj) = New_Data(1);
            CW_T(rowIndex,jj) = New_Data(2);
            Error_Scope(rowIndex,jj) = New_Data(3);
            KddError_N13(rowIndex,jj) = New_Data(4);
            KpError_N10(rowIndex,jj) = New_Data(5);
            Error_N3(rowIndex,jj) = New_Data(6);
            Kp_scaler(rowIndex,jj) = New_Data(7);
            Kd_scaler(rowIndex,jj) = New_Data(8);
            Kp_syn_gain(rowIndex,jj) = New_Data(9);
            Kd_syn_gain(rowIndex,jj) = New_Data(10);
            kd_pre_gain(rowIndex,jj) = New_Data(11);
            
                NNoutput(rowIndex,jj) = ((CCW_T(rowIndex,jj)-Veq)-(CW_T(rowIndex,jj)-Veq))*convert_to_nA; % total control torque from AM, adjusted by equilibrium voltage multiplied by amplifier gain

%                             if abs(NNoutput(rowIndex,jj) - NNoutput(rowIndex-1,jj)) > 5
%                                 fprintf('Serial Error in Frequency %d, rowIndex %d \n',jj,rowIndex)
%                             end

                if rowIndex > delaySteps
                    Tapp(rowIndex) = NNoutput(rowIndex-delaySteps,jj)*gain; %convert theta_out (input to controller to current
                end
%                         Tapp(rowIndex) = NNoutput(rowIndex,jj)*gain;
%                         NNinput_des(rowIndex,jj) = amp(bb)*sin(Input_Frequency*time(rowIndex,jj));  % input amp in nA; input freq in rad/s

                    % Inverted Pendulum Simulation
                    [theta(rowIndex,:)] = Simulate_InvertedPend_Tiff(theta(rowIndex-1,:),Tapp(rowIndex),dt,dt,1); %[position, velocity] inputs(current position/velocity, current corrective Torque applied, simulation time, dt, howmuch) for howmuch = 1, simulation returns only final calculated value, so with time = dt, simulation only runs through one time each time it's accessed in the main code
                    NNinput_angle(rowIndex,jj) = theta(rowIndex,1); 

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
        fprintf('simulation done \n');
        
        Error_N3(:,jj) = (Error_N3(:,jj) - Veq)*1000;
        Error_Scope(:,jj) = (Error_Scope(:,jj) - Veq)*1000;
        KpError_N10(:,jj) = (KpError_N10(:,jj) - Veq)*1000;
        KddError_N13(:,jj) = (KddError_N13(:,jj) - Veq)*1000;
        Kp_scaler(:,jj) = (Kp_scaler(:,jj) - Veq)*1000;
        Kd_scaler(:,jj) = (Kd_scaler(:,jj) - Veq)*1000;
        Kp_syn_gain(:,jj) = (Kp_syn_gain(:,jj) - Veq)*1000;
        Kd_syn_gain(:,jj) = (Kd_syn_gain(:,jj) - Veq)*1000;
        kd_pre_gain(:,jj) = (kd_pre_gain(:,jj) - Veq)*1000;
    
    end
    
    
    
    timeNNinput =       [time(1:ntime-1)    NNinput_des(1:ntime-1)   ];
    timeError_N3 =      [time(1:ntime-1)    Error_N3(1:ntime-1)      ];
    timeError_Scope =   [time(1:ntime-1)    Error_Scope(1:ntime-1)   ];
    timeKpError_N10 =   [time(1:ntime-1)    KpError_N10(1:ntime-1)   ];
    timeKddError_N13 =  [time(1:ntime-1)    KddError_N13(1:ntime-1)  ];
    timeNNoutput =      [time(1:ntime-1)    NNoutput(1:ntime-1)      ];
    timeNNangleOut =    [time(1:ntime-1)    NNinput_angle(1:ntime-1) ];  % actual position in nA
    timeKp_scaler =     [time(1:ntime-1)    Kp_scaler(1:ntime-1)     ];
    timeKd_scaler =     [time(1:ntime-1)    Kd_scaler(1:ntime-1)     ];
    timeKp_syn_gain =   [time(1:ntime-1)    Kp_syn_gain(1:ntime-1)   ];
    timeKd_syn_gain =   [time(1:ntime-1)    Kd_syn_gain(1:ntime-1)   ];
    timekd_pre_gain =   [time(1:ntime-1)    kd_pre_gain(1:ntime-1)   ];

    % UPDATE FORMAT SPEC TO #FREQ + 1 AND UPDATE FILENAME
        formatSpec = '%f %f\n'; % for 1 input
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\NNinput_WNCL_50(',num2str(amp(bb)),').dat'],'w');
    [row,c] = size(timeNNinput);
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeNNinput(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\ErrorN3_WNCL_50(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeError_N3(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\ErrorScope_WNCL_50(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeError_Scope(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\Kp_WNCL_50(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeKpError_N10(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\Kd_WNCL_50(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeKddError_N13(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\NNOutput_WNCL_50(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeNNoutput(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\NNangleOut_WNCL_50(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeNNangleOut(zz,:));
    end
    
    fclose(fid(bb));

    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\Kpscaler_WNCL_50(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeKp_scaler(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\Kdscaler_WNCL_50(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeKd_scaler(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\Kpsyngain_WNCL_50(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeKp_syn_gain(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\Kdsyngain_WNCL_50(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeKd_syn_gain(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\kdpregain_WNCL_50(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timekd_pre_gain(zz,:));
    end
    
    fclose(fid(bb));    
    
end

% loopData = statistics(rateObj)
% delete(rateObj)

%% time domain plots for sanity checks

meanError = mean(abs(Error_N3))
meanTorque = mean(abs(NNoutput))
maxError = max(abs(Error_N3))

figure
plot(time(1:ntime-1),Error_N3(1:ntime-1),time(1:ntime-1),Error_Scope(1:ntime-1))
% plot(time(1:ntime-1),Error_Scope(1:ntime-1))
title('Error Signal')
legend('Error N3','Error Scope')
% xlim([0 800])
% ylim([-10 10])

figure
hold on
plot(time(1:ntime-1),KpError_N10(1:ntime-1))
plot(time(1:ntime-1),Kp_syn_gain(1:ntime-1))
title('Kp*Error Signal')
legend('Kp*Error Scope','Kp syn gain')
ylim([-1 10])

figure
hold on
plot(time(1:ntime-1),KddError_N13(1:ntime-1))
plot(time(1:ntime-1),kd_pre_gain(1:ntime-1))
plot(time(1:ntime-1),Kd_syn_gain(1:ntime-1))
title('Kd*dError Signal')
legend('Kd*de/dt Scope - syn gain & mod','kd pre gain','Kd syn gain only')
ylim([-1 25])

figure
plot(time(1:ntime-1),NNinput_des(1:ntime-1))
hold on
plot(time(1:ntime-1),NNoutput(1:ntime-1))
title('PD Input Output')
legend('NNinput','NNoutput')
ylim([-25 25])

figure
plot(time(1:ntime-1),NNinput_angle(1:ntime-1))
title('System Output')
ylim([-20 20])

figure
hold on
plot(time(1:ntime-1),Kp_scaler(1:ntime-1),'o')
plot(time(1:ntime-1),Kd_scaler(1:ntime-1),'x')
legend('Kp scaler','Kd scaler')

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