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
delay = .1;                % s
delaySteps = delay/dt;      % delay must be divisible into integers by dt
% gain = 1388;  % amplifier gain out of controller to simulation
gain = 833;  % amplifier gain out of controller to simulation
% rad_to_nA = 1.5*(6/pi()); % map radians to current, 30 deg = pi/6 = 1.5 nA

loopTime = zeros(n,1);

% setup serial package for writing to animatlab
pos_head = [255 255 1 18 0 23 0]; %defines header, message ID (always 1), message size (2 bytes) and Data ID for Actual Angle = 23 in AM (2 bytes), starting actual angle = 0
pos_des_head = [24 0]; %[data ID for desired angle = 24 in Animatlab Input Serial IO properties, desired position = 0]
pos_msg = zeros(1,12);

% auto bode frequencies and input amplitudes to loop through
% amp = [0.75 1 1.5 2 3 4];
amp = [1 2];
% amp = [2];  % input amplitude nA
N_amps = length(amp);
convert_to_nA = 1000;
% frequencies_Hz = [1 3 5 10 20 25 40 50 80 100];    % if getting errors, check that Ts divisible by all freqs for auto bode
% frequencies_Hz = [0.05 0.1 0.16 0.2 0.25 0.32 0.5 1 2 2.5];       % 10
% frequencies_Hz = [0.25 0.6 1 1.25 1.5 2 2.5 4 10 20];     % 10
% frequencies_Hz = [0.1 0.16 0.2 0.25 0.4 0.5 1 1.6 2];     % 9 frequencies
% frequencies_Hz = [0.16 0.2 0.25 0.4 0.5 1 1.6 2];     % 8 frequencies
% frequencies_Hz = [0.16 0.2 0.25 0.4 0.5 1 1.6 2 2.25];     % 9 frequencies
% frequencies_Hz = [0.16 0.2 0.25 0.32 0.4 0.5 0.65 0.8 1 1.6 2 2.25];     % 12 frequencies
frequencies_Hz = [1 2 3];
frequencies = 2*pi*frequencies_Hz;

% % set up loop timer (robotics toolbox)
% desiredRate = Ts;
% rateObj = robotics.Rate(desiredRate);
% rateObj.OverrunAction  = 'slip';

% reset(rateObj);

% Auto Bode settings
N_frequencies      = length(frequencies);
N_test_cycles                = 15;
% N_points_per_cycle           = 30;

N_mean_percent               = 20;  % percent of highest and lowest amps and phases from a single cycle which are disregarded as outliers
N_initial_discard_percent    = 30;  % percent of initial cycles that are tossed out initially
    
%% Pull Torques from AM, run IP simulation, Send Actual Position to AM

for bb=1:N_amps

    phase_save_p = zeros(1,N_frequencies);
    Amplitude_save_p = zeros(1,N_frequencies);
    Amplitude_ratio_stdev_p = zeros(1,N_frequencies);
    
    phase_save_p_sys = zeros(1,N_frequencies);
    Amplitude_save_p_sys = zeros(1,N_frequencies);
    Amplitude_ratio_stdev_p_sys = zeros(1,N_frequencies);
    
    phase_save_p_Kd = zeros(1,N_frequencies);
    Amplitude_save_ratio_Kd = zeros(1,N_frequencies);
    Amplitude_ratio_stdev_p_Kd = zeros(1,N_frequencies);

    phase_save_p_Kp = zeros(1,N_frequencies);
    Amplitude_save_ratio_Kp = zeros(1,N_frequencies);
    Amplitude_ratio_stdev_p_Kp = zeros(1,N_frequencies);
    
    Amplitude_Out_save_Kd = zeros(1,N_frequencies);
    Amplitude_in_save_Kd = zeros(1,N_frequencies);
    Amplitude_out_save_Kd_stdev = zeros(1,N_frequencies);
    Amplitude_in_save_Kd_stdev = zeros(1,N_frequencies);
    
    Amplitude_Out_save_Kp = zeros(1,N_frequencies);
    Amplitude_in_save_Kp = zeros(1,N_frequencies);
    Amplitude_out_save_Kp_stdev = zeros(1,N_frequencies);
    Amplitude_in_save_Kp_stdev = zeros(1,N_frequencies);
    
    phase_save_p_Kpsyn = zeros(1,N_frequencies);
    Amplitude_save_p_Kpsyn = zeros(1,N_frequencies);
    Amplitude_ratio_stdev_p_Kpsyn = zeros(1,N_frequencies);
    
    phase_save_p_Kdsyn = zeros(1,N_frequencies);
    Amplitude_save_p_Kdsyn = zeros(1,N_frequencies);
    Amplitude_ratio_stdev_p_Kdsyn = zeros(1,N_frequencies);
    
    phase_save_p_kdbase = zeros(1,N_frequencies);
    Amplitude_save_p_kdbase = zeros(1,N_frequencies);
    Amplitude_ratio_stdev_p_kdbase = zeros(1,N_frequencies);
    
     
    
%     NNoutput = zeros(N_test_cycles*N_points_per_cycle,N_frequencies);    % sum of NN output
%     NNinput_des = zeros(N_test_cycles*N_points_per_cycle,N_frequencies); % 
%     time = zeros(N_test_cycles*N_points_per_cycle+1,N_frequencies);

    NNoutput = zeros(n,N_frequencies);    % sum of NN output
    NNinput_des = zeros(n,N_frequencies); %
    NNinput_angle = zeros(n,N_frequencies); % position as current to apply to AM (zero to N1)
%     NNinput_angle(1,:) = 0.1*rad_to_nA;     % initial displacement ~5.7 degrees
    time = zeros(n,N_frequencies);
    CW_T = zeros(n,N_frequencies);          % CW torque commmand vector
    CCW_T = zeros(n,N_frequencies);         % CCW torque command vector
    Error_Scope = zeros(n,N_frequencies);      % Error Neuron 3 membrane voltage
    KpError_N10 = zeros(n,N_frequencies);   % Kp*Error Neuron 10 membrane voltage
    KddError_N13 = zeros(n,N_frequencies);  % Kp*dError Neuron 13 membrane voltage
    Error_N3 = zeros(n,N_frequencies);      % Error Neuron 3 membrane voltage
    Kp_scaler = zeros(n,N_frequencies);      % Kp scaler N8 membrane voltage
    Kd_scaler = zeros(n,N_frequencies);      % Kd scaler N11 membrane voltage
    Kp_syn_gain = zeros(n,N_frequencies);      % Kp syn gain membrane voltage
    Kd_syn_gain = zeros(n,N_frequencies);      % Kd syn gain membrane voltage
    kd_pre_gain = zeros(n,N_frequencies);      % kd*de/dt CW membrane voltage
    
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

        theta = zeros(n,2); % position (radians) and velocity (rad/s)
        Tapp = zeros(n,1);  % Torque applied (controller output*gain)
        Input_Frequency     = frequencies(jj);          %Rad/sec
        P                   = 2*pi/Input_Frequency;     %Period (sec)
        Input_Frequency_Hz  = frequencies_Hz(jj);
        N_points_per_cycle  = floor(Sample_Frequency/Input_Frequency_Hz);
        N_points_per_freq   = N_points_per_cycle*N_test_cycles;
%         dt                  = P/N_points_per_cycle;     %Sample time (sec)
        tend                = N_test_cycles*P;          %End of cycle time (sec)
        rowIndex = 2;
        
        fprintf('start frequency %i Hz \n',frequencies_Hz(jj));       

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

                    if abs(NNoutput(rowIndex,jj) - NNoutput(rowIndex-1,jj)) > 5
                        fprintf('Serial Error in Frequency %d, rowIndex %d \n',jj,rowIndex)
                    end

                if rowIndex > delaySteps
                    Tapp(rowIndex) = NNoutput(rowIndex-delaySteps,jj)*gain; %convert theta_out (input to controller to current
                end
%                         Tapp(rowIndex) = NNoutput(rowIndex,jj)*gain;
                NNinput_des(rowIndex,jj) = amp(bb)*sin(Input_Frequency*time(rowIndex,jj));  % input amp in nA; input freq in rad/s

                    % Inverted Pendulum Simulation
                    [theta(rowIndex,:)] = Simulate_InvertedPend_Tiff(theta(rowIndex-1,:),Tapp(rowIndex),dt,dt,1); %[position, velocity] inputs(current position/velocity, current corrective Torque applied, simulation time, dt, howmuch) for howmuch = 1, simulation returns only final calculated value, so with time = dt, simulation only runs through one time each time it's accessed in the main code
                    NNinput_angle(rowIndex,jj) = theta(rowIndex,1); %convert theta_out (input to controller to current

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

%         m = N_points_per_cycle*N_test_cycles + 1;
        
        Error_N3(:,jj) = (Error_N3(:,jj) - Veq)*1000;
        Error_Scope(:,jj) = (Error_Scope(:,jj) - Veq)*1000;
        KpError_N10(:,jj) = (KpError_N10(:,jj) - Veq)*1000;
        KddError_N13(:,jj) = (KddError_N13(:,jj) - Veq)*1000;
        Kp_scaler(:,jj) = (Kp_scaler(:,jj) - Veq)*1000;
        Kd_scaler(:,jj) = (Kd_scaler(:,jj) - Veq)*1000;
        Kp_syn_gain(:,jj) = (Kp_syn_gain(:,jj) - Veq)*1000;
        Kd_syn_gain(:,jj) = (Kd_syn_gain(:,jj) - Veq)*1000;
        kd_pre_gain(:,jj) = (kd_pre_gain(:,jj) - Veq)*1000;
        
        % Clear all Ouputs and Auto Bode
        % open loop PD bode
        clear Amplitude_save phase_save

        [Amplitude_save,phase_save,Amplitude_ratio_stdev] = Auto_Bode(Error_N3(1:N_points_per_freq,jj),NNoutput(1:N_points_per_freq,jj),N_test_cycles,N_points_per_cycle,N_mean_percent,N_initial_discard_percent);

        phase_save_p(jj)     = phase_save;
        Amplitude_save_p(jj) = Amplitude_save;
        Amplitude_ratio_stdev_p(jj) = Amplitude_ratio_stdev;
        
        % closed loop system bode
        clear Amplitude_save_sys phase_save_sys
        [Amplitude_save_sys,phase_save_sys,Amplitude_ratio_stdev_sys] = Auto_Bode(NNinput_des(1:N_points_per_freq,jj),NNinput_angle(1:N_points_per_freq,jj),N_test_cycles,N_points_per_cycle,N_mean_percent,N_initial_discard_percent);

        phase_save_p_sys(jj)     = phase_save_sys;
        Amplitude_save_p_sys(jj) = Amplitude_save_sys;
        Amplitude_ratio_stdev_p_sys(jj) = Amplitude_ratio_stdev_sys;
        
        % Kd Auto Bode
        clear Amplitude_ratio_Kd phase_save_Kd Amplitude_out_Kd Amplitude_in_Kd
        
        [Amplitude_ratio_Kd,Amplitude_out_Kd,Amplitude_out_Kd_stdev,Amplitude_in_Kd,Amplitude_in_Kd_stdev,phase_save_Kd,Amplitude_ratio_stdev_Kd] = Auto_Bode_Tiff(Error_N3(1:N_points_per_freq,jj),2*KddError_N13(1:N_points_per_freq,jj),N_test_cycles,N_points_per_cycle,N_mean_percent,N_initial_discard_percent);

        phase_save_p_Kd(jj)     = phase_save_Kd;
        Amplitude_save_ratio_Kd(jj) = Amplitude_ratio_Kd;
        Amplitude_ratio_stdev_p_Kd(jj) = Amplitude_ratio_stdev_Kd;
        Amplitude_Out_save_Kd(jj) = Amplitude_out_Kd/2; % divide by 2 to convert from peak to peak amplitudes to standard amplitudes
        Amplitude_out_save_Kd_stdev(jj) = Amplitude_out_Kd_stdev/2;
        Amplitude_in_save_Kd(jj) = Amplitude_in_Kd/2;
        Amplitude_in_save_Kd_stdev(jj) = Amplitude_in_Kd_stdev/2;

        % Kp Auto Bode
        clear Amplitude_ratio_Kp phase_save_Kp Amplitude_out_Kp Amplitude_in_Kp
        
        [Amplitude_ratio_Kp,Amplitude_out_Kp,Amplitude_out_Kp_stdev,Amplitude_in_Kp,Amplitude_in_Kp_stdev,phase_save_Kp,Amplitude_ratio_stdev_Kp] = Auto_Bode_Tiff(Error_N3(1:N_points_per_freq,jj),2*KpError_N10(1:N_points_per_freq,jj),N_test_cycles,N_points_per_cycle,N_mean_percent,N_initial_discard_percent);

        phase_save_p_Kp(jj)     = phase_save_Kp;
        Amplitude_save_ratio_Kp(jj) = Amplitude_ratio_Kp;
        Amplitude_ratio_stdev_p_Kp(jj) = Amplitude_ratio_stdev_Kp;
        Amplitude_Out_save_Kp(jj) = Amplitude_out_Kp/2;
        Amplitude_out_save_Kp_stdev(jj) = Amplitude_out_Kp_stdev/2;
        Amplitude_in_save_Kp(jj) = Amplitude_in_Kp/2;
        Amplitude_in_save_Kp_stdev(jj) = Amplitude_in_Kp_stdev/2;
        
        % Kp synaptic gain/Error_N3
        clear Amplitude_save_Kpsyn phase_save_Kpsyn

        [Amplitude_save_Kpsyn,phase_save_Kpsyn,Amplitude_ratio_stdev_Kpsyn] = Auto_Bode(Error_N3(1:N_points_per_freq,jj),2*Kp_syn_gain(1:N_points_per_freq,jj),N_test_cycles,N_points_per_cycle,N_mean_percent,N_initial_discard_percent);
            % the above is more accurate than Kpsyn/ErrorScope, because
            % that divides out the base gain of 1.11 and doesn't reflect
            % the effect of varied inputs
        phase_save_p_Kpsyn(jj)     = phase_save_Kpsyn;
        Amplitude_save_p_Kpsyn(jj) = Amplitude_save_Kpsyn;
        Amplitude_ratio_stdev_p_Kpsyn(jj) = Amplitude_ratio_stdev_Kpsyn;
        
        % Kd synaptic gain/kd pre gain
        clear Amplitude_save_Kdsyn phase_save_Kdsyn

        [Amplitude_save_Kdsyn,phase_save_Kdsyn,Amplitude_ratio_stdev_Kdsyn] = Auto_Bode(kd_pre_gain(1:N_points_per_freq,jj),Kd_syn_gain(1:N_points_per_freq,jj),N_test_cycles,N_points_per_cycle,N_mean_percent,N_initial_discard_percent);

        phase_save_p_Kdsyn(jj)     = phase_save_Kdsyn;
        Amplitude_save_p_Kdsyn(jj) = Amplitude_save_Kdsyn;
        Amplitude_ratio_stdev_p_Kdsyn(jj) = Amplitude_ratio_stdev_Kdsyn;
        
        % Baseline kd gain = kd pre gain/Error Scope
        clear Amplitude_save_Kdsyn phase_save_Kdsyn

        [Amplitude_save_kdbase,phase_save_kdbase,Amplitude_ratio_stdev_kdbase] = Auto_Bode(Error_Scope(1:N_points_per_freq,jj),kd_pre_gain(1:N_points_per_freq,jj),N_test_cycles,N_points_per_cycle,N_mean_percent,N_initial_discard_percent);

        phase_save_p_kdbase(jj)     = phase_save_kdbase;
        Amplitude_save_p_kdbase(jj) = Amplitude_save_kdbase;
        Amplitude_ratio_stdev_p_kdbase(jj) = Amplitude_ratio_stdev_kdbase;

    end
    
    fprintf('amplitude %d done \n',amp(bb));
    
%     frequencies = frequencies/(2*pi);
%     Amplitude_save_p = 20*log10(Amplitude_save_p);

    % Open Loop PD Bode
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\PD_ABCL_13(',num2str(amp(bb)),').dat'],'w');

    for ii=1:length(frequencies)
        fprintf(fid(bb),'%12e   %12e   %12e   %12e\n' ,frequencies(ii),Amplitude_save_p(ii),phase_save_p(ii),Amplitude_ratio_stdev_p(ii));
    end
    
    % Closed Loop System Bode
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\SYS_ABCL_13(',num2str(amp(bb)),').dat'],'w');

    for ii=1:length(frequencies)
        fprintf(fid(bb),'%12e   %12e   %12e   %12e\n' ,frequencies(ii),Amplitude_save_p_sys(ii),phase_save_p_sys(ii),Amplitude_ratio_stdev_p_sys(ii));
    end
    
    fclose(fid(bb));

    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kd_ABCL_13(',num2str(amp(bb)),').dat'],'w');

    for ii=1:length(frequencies)
        fprintf(fid(bb),'%12e   %12e   %12e   %12e   %12e   %12e   %12e   %12e\n' ,frequencies(ii),Amplitude_save_ratio_Kd(ii),phase_save_p_Kd(ii),Amplitude_in_save_Kd(ii),Amplitude_in_save_Kd_stdev(ii),Amplitude_Out_save_Kd(ii),Amplitude_out_save_Kd_stdev(ii),Amplitude_ratio_stdev_p_Kd(ii));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kp_ABCL_13(',num2str(amp(bb)),').dat'],'w');

    for ii=1:length(frequencies)
        fprintf(fid(bb),'%12e   %12e   %12e   %12e   %12e   %12e   %12e   %12e\n' ,frequencies(ii),Amplitude_save_ratio_Kp(ii),phase_save_p_Kp(ii),Amplitude_in_save_Kp(ii),Amplitude_in_save_Kp_stdev(ii),Amplitude_Out_save_Kp(ii),Amplitude_out_save_Kp_stdev(ii),Amplitude_ratio_stdev_p_Kp(ii));
    end  
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kpsyngain_ABCL_13(',num2str(amp(bb)),').dat'],'w');

    for ii=1:length(frequencies)
        fprintf(fid(bb),'%12e   %12e   %12e   %12e\n' ,frequencies(ii),Amplitude_save_p_Kpsyn(ii),phase_save_p_Kpsyn(ii),Amplitude_ratio_stdev_p_Kpsyn(ii));
    end  
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kdsyngain_ABCL_13(',num2str(amp(bb)),').dat'],'w');

    for ii=1:length(frequencies)
        fprintf(fid(bb),'%12e   %12e   %12e   %12e\n' ,frequencies(ii),Amplitude_save_p_Kdsyn(ii),phase_save_p_Kdsyn(ii),Amplitude_ratio_stdev_p_Kdsyn(ii));
    end  
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\kdbase_ABCL_13(',num2str(amp(bb)),').dat'],'w');

    for ii=1:length(frequencies)
        fprintf(fid(bb),'%12e   %12e   %12e   %12e\n' ,frequencies(ii),Amplitude_save_p_kdbase(ii),phase_save_p_kdbase(ii),Amplitude_ratio_stdev_p_kdbase(ii));
    end  
    
    fclose(fid(bb));
    
    timeNNinput =       [time(1:50000,1)    NNinput_des(1:50000,:)  ];  % desired position input in nA
    timeError_N3 =      [time(1:50000,1)    Error_N3(1:50000,:)     ];  % units membrane voltage (mV) (or nA current), normalized to equ voltage
    timeError_Scope =   [time(1:50000,1)    Error_Scope(1:50000,:)  ];
    timeKpError_N10 =   [time(1:50000,1)    KpError_N10(1:50000,:)  ];  % units membrane voltage (mV) (or nA current), normalized to equ voltage
    timeKddError_N13 =  [time(1:50000,1)    KddError_N13(1:50000,:) ];  % units membrane voltage (mV) (or nA current), normalized to equ voltage
    timeNNoutput =      [time(1:50000,1)    NNoutput(1:50000,:)     ];  % output torque in nA
    timeNNangleOut =    [time(1:50000,1)    NNinput_angle(1:50000,:)];  % actual position in nA
    timeKp_scaler =     [time(1:50000,1)    Kp_scaler(1:50000,:)    ];
    timeKd_scaler =     [time(1:50000,1)    Kd_scaler(1:50000,:)    ];
    timeKp_syn_gain =   [time(1:50000,1)    Kp_syn_gain(1:50000,:)  ];
    timeKd_syn_gain =   [time(1:50000,1)    Kd_syn_gain(1:50000,:)  ];
    timekd_pre_gain =   [time(1:50000,1)    kd_pre_gain(1:50000,:)  ];

    % UPDATE FORMAT SPEC TO #FREQ + 1 AND UPDATE FILENAME
%         formatSpec = '%f %f %f %f %f %f %f %f %f %f %f\n';     % for 10 input frequencies
%         formatSpec = '%f %f %f\n'; % for 2 input frequencies
%         formatSpec = '%f %f %f %f %f\n'; % for 4 input frequencies
%         formatSpec = '%f %f %f %f %f %f %f\n'; % for 6 input frequencies         
%         formatSpec = '%f %f %f %f %f %f %f %f\n'; % for 7 input frequencies
%         formatSpec = '%f %f %f %f %f %f %f %f %f\n'; % for 8 input frequencies
%         formatSpec = '%f %f %f %f %f %f %f %f %f %f\n'; % for 9 input frequencies
%         formatSpec = '%f %f %f %f %f %f %f %f %f %f %f\n'; % for 10 input frequencies
        formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f %f\n'; % for 12 input frequencies
%         formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n'; % for 14 input frequencies
%         formatSpec = '%f %f\n'; % for 1 input frequency
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_NNOutputs\NNinput_ABCL_13(',num2str(amp(bb)),').dat'],'w');
    [row,c] = size(timeNNinput);
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeNNinput(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_NNOutputs\ErrorN3_ABCL_13(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeError_N3(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_NNOutputs\ErrorScope_ABCL_13(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeError_Scope(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_NNOutputs\Kp_ABCL_13(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeKpError_N10(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_NNOutputs\Kd_ABCL_13(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeKddError_N13(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_NNOutputs\NNOutput_ABCL_13(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeNNoutput(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_NNOutputs\NNangleOut_ABCL_13(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeNNangleOut(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_NNOutputs\Kpscaler_ABCL_13(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeKp_scaler(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_NNOutputs\Kdscaler_ABCL_13(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeKd_scaler(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_NNOutputs\Kpsyngain_ABCL_13(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeKp_syn_gain(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_NNOutputs\Kdsyngain_ABCL_13(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeKd_syn_gain(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_NNOutputs\kdpregain_ABCL_13(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timekd_pre_gain(zz,:));
    end
    
    fclose(fid(bb));
    
end

fprintf('simulation done \n');

% loopData = statistics(rateObj)
% delete(rateObj)

%% time domain plots for sanity checks

figure
plot(time(:,1),Error_N3(:,1),time(:,1),Error_Scope(:,1))
title('Error Signal')
legend('Error N3','Error Scope')
ylim([-10 10])

figure
hold on
plot(time(:,1),KpError_N10(:,1))
plot(time(:,1),Kp_syn_gain(:,1))
title('Kp*Error Signal')
legend('Kp*Error Scope','Kp syn gain')
ylim([-2 2])

figure
hold on
plot(time(:,1),KddError_N13(:,1))
plot(time(:,1),kd_pre_gain(:,1))
plot(time(:,1),Kd_syn_gain(:,1))
title('Kd*dError Signal')
legend('Kd*de/dt Scope - syn gain & mod','kd pre gain','Kd syn gain only')
ylim([-3 3])

figure
plot(time(:,1),NNoutput(:,1))
title('NN Output Signal')
ylim([-10 10])
%%
figure
hold on
plot(time(:,1),NNinput_des(:,1))
plot(time(:,1),NNinput_angle(:,1))
plot(time(:,1),NNoutput(:,1))
legend('Command Angle','Actual Angle','NNoutput')
ylim([-2 2])

figure
hold on
plot(time(:,1),Tapp(:,1))
legend('Torque Applied')
ylim([-2000 2000])
%%
figure
hold on
plot(time(:,1),Kp_scaler(:,1),'o')
plot(time(:,1),Kd_scaler(:,1),'x')
legend('Kp scaler','Kd scaler')
% ylim([-2 2])

%% STILL OPEN LOOP PLOTS HERE
s = tf('s');

% Peterka's "normal subject" constants
m = 83.3; %kg 
h = .896; %m - center of mass height
I = 81.1; %kg m^2
g = 9.81;% m/s

figure
semilogx(frequencies_Hz,Amplitude_save_p_sys,'-o',frequencies_Hz,(Amplitude_save_p_sys + 2*Amplitude_ratio_stdev_p_sys),'--',frequencies_Hz,(Amplitude_save_p_sys - 2*Amplitude_ratio_stdev_p_sys),'--')
title('Amplitude')
xlabel('Frequencies Hz')
ylabel('Mag Ratio (abs)')
legend('Magnitude','+ 2std','- 2std','Location','northwest')

figure 
semilogx(frequencies_Hz,phase_save_p_sys,'-o')
title('Phase')
xlabel('Frequencies Hz')
ylabel('Phase degrees')
% legend('Actual','Theoretical','Location','northeast')
% ylim([-100 10])

figure
semilogx(frequencies_Hz,Amplitude_save_ratio_Kp,'m-+',frequencies_Hz,(Amplitude_save_ratio_Kp + 2*Amplitude_ratio_stdev_p_Kp),'--',frequencies_Hz,(Amplitude_save_ratio_Kp - 2*Amplitude_ratio_stdev_p_Kp),'--')
title('Kp Out/Error In')
xlabel('Frequencies Hz')
ylabel('Mag Ratio (abs)')

figure
semilogx(frequencies_Hz,Amplitude_save_ratio_Kd,'k-d',frequencies_Hz,(Amplitude_save_ratio_Kd + 2*Amplitude_ratio_stdev_p_Kd),'--',frequencies_Hz,(Amplitude_save_ratio_Kd - 2*Amplitude_ratio_stdev_p_Kd),'--')
title('Kd*dE Out/Error In')
xlabel('Frequencies Hz')
ylabel('Mag Ratio (dB)')

figure
semilogx(frequencies_Hz,Amplitude_save_p_Kpsyn,'-o',frequencies_Hz,(Amplitude_save_p_Kpsyn + 2*Amplitude_ratio_stdev_p_Kpsyn),'--',frequencies_Hz,(Amplitude_save_p_Kpsyn - 2*Amplitude_ratio_stdev_p_Kpsyn),'--')
hold on
semilogx(frequencies_Hz,Amplitude_save_p_Kdsyn,'-x',frequencies_Hz,(Amplitude_save_p_Kdsyn + 2*Amplitude_ratio_stdev_p_Kdsyn),'--',frequencies_Hz,(Amplitude_save_p_Kdsyn - 2*Amplitude_ratio_stdev_p_Kdsyn),'--')
title('Amplitude')
xlabel('Frequencies Hz')
ylabel('Mag Ratio (abs)')
legend('Kp syn gain','Kd syn gain','Location','northwest')

figure
semilogx(frequencies_Hz,Amplitude_save_p_kdbase,'-o',frequencies_Hz,(Amplitude_save_p_kdbase + 2*Amplitude_ratio_stdev_p_kdbase),'--',frequencies_Hz,(Amplitude_save_p_kdbase - 2*Amplitude_ratio_stdev_p_kdbase),'--')
title('Amplitude')
xlabel('Frequencies Hz')
ylabel('kd Base Gain kd*de/dt CW/Error Scope')
legend('Kp syn gain','Kd syn gain','Location','northwest')

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