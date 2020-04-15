%% Tiffany Hamstreet - Portland State University
%% Animatlab Serial interface with MATLAB (credit to Wade original HEBI/Matlab/AM code)
%% For automatic bode creation with OPEN LOOP PD Controller, includes serial connections to Error (N3), Kp Scope, and Kd Scope Neurons

clear all;
clc;

% setup serial connection properties
delete(instrfindall)
r = serial('COM2');     % receive data from AM over COM1
set(r,'BaudRate',256000); % it is critical that strict baudrate emulation is enabled on the virtual serial port driver
r.Timeout = 20; % default 10
fopen(r);
s = serial('COM3');     % send data to AM over COM2
set(s,'BaudRate',256000); % Animatlab seems to be tolerant of very high virtual baudrates. I could have gone to 10^6
fopen(s);

% initialize variables for serial communication
n = 200000; % initialize the maximum number of serial bytes that will be read before vectors are full (set to well above what sim time needs)

% initialize simulation variables and parameters
Veq = -.060; %V used to normalize commanded torque

Sample_Frequency = 400;  % sample frequency Hz
dt = 1/Sample_Frequency; %s
loopTime = zeros(n,1);
gain = 1;  % amplifier gain out of controller to simulation
pos_head = [255 255 1 18 0 23 0]; %defines header, message ID (always 1), message size (2 bytes) and Data ID for Actual Angle = 23 in AM (2 bytes), starting actual angle = 0
pos_des_head = [24 0]; %[data ID for desired angle = 24 in Animatlab Input Serial IO properties, desired position = 0]
pos_msg = zeros(1,12);
NNinput_angle = zeros(1,n); % position as current to apply to AM (zero to N1)

% rad_to_nA = 1.5*(6/pi());

% auto bode frequencies and input amplitudes to loop through
amp = [0.1 0.8 4];
% amp = [0.1 0.25 0.5 0.75 1 1.5 2 3 4];    % nA
N_amps = length(amp);
convert_to_nA = 1000;
% frequencies = 2*pi*logspace(0,2,15);   %rad/s
% frequencies = [2*pi*100];   %rad/s
% frequencies_Hz = [1 3 5 10 20 25 40 50 80 100];    % if getting errors, check that Ts divisible by all freqs for auto bode
frequencies_Hz = [0.5 1 2 2.5 3 5 10 25 50 80 100];   % 11 frequencies
% frequencies_Hz = [1 2];
frequencies = 2*pi*frequencies_Hz;

% Auto Bode settings
N_frequencies      = length(frequencies);
N_test_cycles                = 15;
% N_points_per_cycle           = 30;

N_mean_percent               = 20;
N_initial_discard_percent    = 20;
    
%% Pull Torques from AM, run IP simulation, Send Actual Position to AM

for bb = 1:N_amps

    phase_save_p = zeros(1,N_frequencies);
    Amplitude_save_p = zeros(1,N_frequencies);
    Amplitude_ratio_stdev_p = zeros(1,N_frequencies);
    
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
    
    for jj = 1:N_frequencies

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
                        
                NNoutput(rowIndex,jj) = ((CCW_T(rowIndex,jj)-Veq)-(CW_T(rowIndex,jj)-Veq))*convert_to_nA*gain; % total control torque from AM, adjusted by equilibrium voltage (subtracts 1 because the cmd_ccw and cmd_cw indexes were already advanced) multiplied by overall gain

%                 % Inverted Pendulum Simulation
%                 [theta(cmd+1,:)] = Simulate_InvertedPend_Tiff(theta(cmd,:),Tapp(cmd),dt,dt,1); %[position, velocity] inputs(current position/velocity, current corrective Torque applied, simulation time, dt, howmuch) for howmuch = 1, simulation returns only final calculated value, so with time = dt, simulation only runs through one time each time it's accessed in the main code
%                 pos_nA(cmd+1) = theta(cmd+1,1)*rad_to_nA; %convert theta_out to current

                NNinput_des(rowIndex,jj) = amp(bb)*sin(Input_Frequency*time(rowIndex,jj));  %input freq in rad/s

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
        
        % Kd Auto Bode
        clear Amplitude_save_Kd phase_save_Kd
        
        [Amplitude_ratio_Kd,Amplitude_out_Kd,Amplitude_out_Kd_stdev,Amplitude_in_Kd,Amplitude_in_Kd_stdev,phase_save_Kd,Amplitude_ratio_stdev_Kd] = Auto_Bode_Tiff(Error_N3(1:N_points_per_freq,jj),2*KddError_N13(1:N_points_per_freq,jj),N_test_cycles,N_points_per_cycle,N_mean_percent,N_initial_discard_percent);

        phase_save_p_Kd(jj)     = phase_save_Kd;
        Amplitude_save_ratio_Kd(jj) = Amplitude_ratio_Kd;
        Amplitude_ratio_stdev_p_Kd(jj) = Amplitude_ratio_stdev_Kd;
        Amplitude_Out_save_Kd(jj) = Amplitude_out_Kd/2; % divide by 2 to convert from peak to peak amplitudes to standard amplitudes
        Amplitude_out_save_Kd_stdev(jj) = Amplitude_out_Kd_stdev/2;
        Amplitude_in_save_Kd(jj) = Amplitude_in_Kd/2;
        Amplitude_in_save_Kd_stdev(jj) = Amplitude_in_Kd_stdev/2;

        % Kp Auto Bode
        clear Amplitude_save_Kp phase_save_Kp
        
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
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\PD_ABOL_12(',num2str(amp(bb)),').dat'],'w');

    for ii=1:length(frequencies)
        fprintf(fid(bb),'%12e   %12e   %12e   %12e\n' ,frequencies(ii),Amplitude_save_p(ii),phase_save_p(ii),Amplitude_ratio_stdev_p(ii));
    end
    
    fclose(fid(bb));

    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kd_ABOL_12(',num2str(amp(bb)),').dat'],'w');

    for ii=1:length(frequencies)
        fprintf(fid(bb),'%12e   %12e   %12e   %12e   %12e   %12e   %12e   %12e\n' ,frequencies(ii),Amplitude_save_ratio_Kd(ii),phase_save_p_Kd(ii),Amplitude_in_save_Kd(ii),Amplitude_in_save_Kd_stdev(ii),Amplitude_Out_save_Kd(ii),Amplitude_out_save_Kd_stdev(ii),Amplitude_ratio_stdev_p_Kd(ii));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kp_ABOL_12(',num2str(amp(bb)),').dat'],'w');

    for ii=1:length(frequencies)
        fprintf(fid(bb),'%12e   %12e   %12e   %12e   %12e   %12e   %12e   %12e\n' ,frequencies(ii),Amplitude_save_ratio_Kp(ii),phase_save_p_Kp(ii),Amplitude_in_save_Kp(ii),Amplitude_in_save_Kp_stdev(ii),Amplitude_Out_save_Kp(ii),Amplitude_out_save_Kp_stdev(ii),Amplitude_ratio_stdev_p_Kp(ii));
    end  
    
    fclose(fid(bb));
    
        fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kpsyngain_ABOL_12(',num2str(amp(bb)),').dat'],'w');

    for ii=1:length(frequencies)
        fprintf(fid(bb),'%12e   %12e   %12e   %12e\n' ,frequencies(ii),Amplitude_save_p_Kpsyn(ii),phase_save_p_Kpsyn(ii),Amplitude_ratio_stdev_p_Kpsyn(ii));
    end  
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kdsyngain_ABOL_12(',num2str(amp(bb)),').dat'],'w');

    for ii=1:length(frequencies)
        fprintf(fid(bb),'%12e   %12e   %12e   %12e\n' ,frequencies(ii),Amplitude_save_p_Kdsyn(ii),phase_save_p_Kdsyn(ii),Amplitude_ratio_stdev_p_Kdsyn(ii));
    end  
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\kdbase_ABOL_12(',num2str(amp(bb)),').dat'],'w');

    for ii=1:length(frequencies)
        fprintf(fid(bb),'%12e   %12e   %12e   %12e\n' ,frequencies(ii),Amplitude_save_p_kdbase(ii),phase_save_p_kdbase(ii),Amplitude_ratio_stdev_p_kdbase(ii));
    end  
    
    fclose(fid(bb));
    
    timeNNinput =       [time(1:50000,1)    NNinput_des(1:50000,:)  ];
    timeError_N3 =      [time(1:50000,1)    Error_N3(1:50000,:)     ];
    timeError_Scope =   [time(1:50000,1)    Error_Scope(1:50000,:)  ];
    timeKpError_N10 =   [time(1:50000,1)    KpError_N10(1:50000,:)  ];  % units membrane voltage (mV) (or nA current), normalized to equ voltage
    timeKddError_N13 =  [time(1:50000,1)    KddError_N13(1:50000,:) ];  % units membrane voltage (mV) (or nA current), normalized to equ voltage
    timeNNoutput =      [time(1:50000,1)    NNoutput(1:50000,:)     ];  % output torque in nA
    timeKp_scaler =     [time(1:50000,1)    Kp_scaler(1:50000,:)    ];
    timeKd_scaler =     [time(1:50000,1)    Kd_scaler(1:50000,:)    ];
    timeKp_syn_gain =   [time(1:50000,1)    Kp_syn_gain(1:50000,:)  ];
    timeKd_syn_gain =   [time(1:50000,1)    Kd_syn_gain(1:50000,:)  ];
    timekd_pre_gain =   [time(1:50000,1)    kd_pre_gain(1:50000,:)  ];
    

    % UPDATE FORMAT SPEC TO #FREQ + 1 AND UPDATE FILENAME
%         formatSpec = '%f %f %f %f %f %f %f %f %f %f %f\n';     % for 10 input frequencies
        formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f\n';     % for 11 input frequencies
%         formatSpec = '%f %f %f %f %f %f %f %f\n';     % for 7 input frequencies
%         formatSpec = '%f %f %f\n'; % for 2 input frequencies
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_NNOutputs\NNinput_ABOL_12(',num2str(amp(bb)),').dat'],'w');
    [row,c] = size(timeNNinput);
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeNNinput(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_NNOutputs\ErrorN3_ABOL_12(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeError_N3(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_NNOutputs\ErrorScope_ABOL_12(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeError_Scope(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_NNOutputs\Kp_ABOL_12(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeKpError_N10(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_NNOutputs\Kd_ABOL_12(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeKddError_N13(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_NNOutputs\NNOutput_ABOL_12(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeNNoutput(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_NNOutputs\Kpscaler_ABOL_12(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeKp_scaler(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_NNOutputs\Kdscaler_ABOL_12(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeKd_scaler(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_NNOutputs\Kpsyngain_ABOL_12(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeKp_syn_gain(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_NNOutputs\Kdsyngain_ABOL_12(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timeKd_syn_gain(zz,:));
    end
    
    fclose(fid(bb));
    
    fid(bb) = fopen(['C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_NNOutputs\kdpregain_ABOL_12(',num2str(amp(bb)),').dat'],'w');
    for zz = 1:row
        fprintf(fid(bb),formatSpec,timekd_pre_gain(zz,:));
    end
    
    fclose(fid(bb));
    
end

fprintf('simulation done \n');

% loopData = statistics(rateObj)
% delete(rateObj)

%% add post-processing for amplitude and phase finding

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

figure
hold on
plot(time(:,1),Kp_scaler(:,1),'o')
plot(time(:,1),Kd_scaler(:,1),'x')
legend('Kp scaler','Kd scaler')
% ylim([-2 2])

%%
s = tf('s');

Kp = 1;
Kd = 0.04;
unknownGain = 1;
tau = 0;

fb1 = 5;         % cutoff frequency in Hz (base 2)
lowFiltDerivative = (1/(s/(2*pi*fb1)+1));  % low pass filter

fb2 = 80;         % cutoff frequency in Hz (base 25)
lowFiltNeuron = (1/(s/(2*pi*fb2)+1));  % low pass filter

[num,den] = pade(tau,1);    % pade delay approximation first order (in order to use impulse command)
delay = tf(num,den);

% sys = unknownGain*(Kp*lowFiltNeuron + Kd*s*lowFiltNeuron*lowFiltDerivative);
sys2 = lowFiltNeuron^4*(Kp + Kd*s*lowFiltDerivative);
sysDelay = delay*sys2;

[magDelay,phaseDelay] = bode(sysDelay,frequencies);    % requires frequency in rad/s
magDelay = squeeze(magDelay);
phaseDelay = squeeze(phaseDelay);

% phase_save_p(1) = phase_save_p(1) + 360;

% [mag,phase] = bode(tf1,frequencies);
% mag = 20*log10(squeeze(mag));
% phase = squeeze(phase);
% 
% opts = bodeoptions('cstprefs');
% opts.FreqUnits = 'Hz';
% figure
% bode(sys*lowFiltDerivative*lowFiltNeuron,opts)
% title('Baseline for Neural Controller: Open Loop PD Bode Kp = 1, Kd = 0.0319')
% xlim([1 100])

figure
semilogx(frequencies_Hz,Amplitude_save_p,'-o',frequencies_Hz,magDelay,'-x')
title('Amplitude')
xlabel('Frequencies Hz')
ylabel('Mag Ratio (abs)')
legend('Actual','Fit','Location','northwest')

figure 
semilogx(frequencies_Hz,phase_save_p,'-o',frequencies_Hz,phaseDelay)
title('Phase')
xlabel('Frequencies Hz')
ylabel('Phase degrees')
% ylim([-100 10])

figure
semilogx(frequencies_Hz,Amplitude_save_ratio_Kp,'m-+')
title('Kp Out/Error In')
xlabel('Frequencies Hz')
ylabel('Mag Ratio (abs)')

figure
semilogx(frequencies_Hz,Amplitude_save_ratio_Kd,'k-d')
hold on
semilogx(frequencies_Hz,Amplitude_save_p_kdbase,'-o')
title('Kd*dE Out/Error In & kd base gain kd*de Out/Error Scope In')
xlabel('Frequencies Hz')
ylabel('Mag Ratio (dB)')
legend('Overall Kd Circuit Gain','kd Baselind Gain')

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

% have to fix time scale to see signal
% figure
% plot(time, NNinput_des)
% xlabel('time (s)')
% ylabel('Input Voltage (mV)')
% title('Input to open loop PD Controller (mV)')
% 
% figure
% plot(time, NNoutput)
% % xlim([0 2])
% xlabel('time (s)')
% ylabel('Output Voltage (mV)')
% title('Output of open loop PD Controller (mV)')