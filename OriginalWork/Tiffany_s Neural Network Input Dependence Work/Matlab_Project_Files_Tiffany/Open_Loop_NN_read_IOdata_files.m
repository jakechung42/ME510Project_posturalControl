% Read data from NNinput and output data files

% MAKE SURE TO UPDATE FILENAMES AND AMPLITUDES IN ()

clear
clc

% Import Data
% assumes all have the same time scale
A = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_NNOutputs\NNinput_ABOL_3(4).dat');   % Desired input to NN
% A = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\NNinputABCLTest1(0.4).dat');   % Desired input to NN
time = A(:,1);      % s
Af1 = A(:,2);  
Af2 = A(:,3); 
Af3 = A(:,4);  
Af4 = A(:,5); 
Af5 = A(:,6);  
Af6 = A(:,7); 
Af7 = A(:,8);  
Af8 = A(:,9); 
Af9 = A(:,10); 
Af10 = A(:,11);  
Af11 = A(:,12); 
% Af12 = A(:,13); 
% Af13 = A(:,14); 

B = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_NNOutputs\ErrorN3_ABOL_3(4).dat');   % Error N3 outputs
% B = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\ErrorN3ABCLTest1(0.4).dat');   % Error N3 outputs
Bf1 = B(:,2);  
Bf2 = B(:,3);
Bf3 = B(:,4);  
Bf4 = B(:,5); 
Bf5 = B(:,6);  
Bf6 = B(:,7); 
Bf7 = B(:,8);  
Bf8 = B(:,9); 
Bf9 = B(:,10); 
Bf10 = B(:,11);  
Bf11 = B(:,12); 
% Bf12 = B(:,13); 
% Bf13 = B(:,14); 

C = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_NNOutputs\Kp_ABOL_3(4).dat');   % Kp*Error circuit outputs
% C = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\KpOutputABCLTest1(0.4).dat');   % Kp*Error circuit outputs
Cf1 = C(:,2);  
Cf2 = C(:,3);
Cf3 = C(:,4);  
Cf4 = C(:,5); 
Cf5 = C(:,6);  
Cf6 = C(:,7); 
Cf7 = C(:,8);  
Cf8 = C(:,9); 
Cf9 = C(:,10); 
Cf10 = C(:,11);  
Cf11 = C(:,12); 
% Cf12 = C(:,13); 
% Cf13 = C(:,14); 

D = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_NNOutputs\Kd_ABOL_3(4).dat');   % Kd*de/dt circuit outputs
% D = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\KdOutputABCLTest1(0.4).dat');   % Kd*de/dt circuit outputs
Df1 = D(:,2);  
Df2 = D(:,3);
Df3 = D(:,4);  
Df4 = D(:,5); 
Df5 = D(:,6);  
Df6 = D(:,7); 
Df7 = D(:,8);  
Df8 = D(:,9); 
Df9 = D(:,10); 
Df10 = D(:,11);  
Df11 = D(:,12); 
% Df12 = D(:,13); 
% Df13 = D(:,14);

E = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_NNOutputs\NNOutput_ABOL_3(4).dat');   % NN circuit outputs
% E = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\NNOutputABCLTest1(0.4).dat');   % NN circuit outputs
Ef1 = E(:,2);  
Ef2 = E(:,3);
Ef3 = E(:,4);  
Ef4 = E(:,5); 
Ef5 = E(:,6);  
Ef6 = E(:,7); 
Ef7 = E(:,8);  
Ef8 = E(:,9); 
Ef9 = E(:,10); 
Ef10 = E(:,11);  
Ef11 = E(:,12); 
% Ef12 = E(:,13); 
% Ef13 = E(:,14); 

G = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_NNOutputs\Kpscaler_ABOL_3(4).dat');   % Kp scaler
% G = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\KpscalerABCLTest1(0.4).dat');   % Kp scaler
Gf1 = G(:,2);  
Gf2 = G(:,3);
Gf3 = G(:,4);  
Gf4 = G(:,5); 
Gf5 = G(:,6);  
Gf6 = G(:,7); 
Gf7 = G(:,8);  
Gf8 = G(:,9); 
Gf9 = G(:,10); 
Gf10 = G(:,11);  
Gf11 = G(:,12); 
% Gf12 = G(:,13); 
% Gf13 = G(:,14); 

H = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_NNOutputs\Kdscaler_ABOL_3(4).dat');   % Kd scaler
% H = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\KdscalerABCLTest1(0.4).dat');   % Kd scaler
Hf1 = H(:,2);  
Hf2 = H(:,3);
Hf3 = H(:,4);  
Hf4 = H(:,5); 
Hf5 = H(:,6);  
Hf6 = H(:,7); 
Hf7 = H(:,8);  
Hf8 = H(:,9); 
Hf9 = H(:,10); 
Hf10 = H(:,11);  
Hf11 = H(:,12); 
% Hf12 = H(:,13); 
% Hf13 = H(:,14); 

J = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_NNOutputs\Kpsyngain_ABOL_3(4).dat');   % Kp syn gain only
% J = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\KpsyngainABCLTest1(0.4).dat');   % Kp syn gain only
Jf1 = J(:,2);  
Jf2 = J(:,3);
Jf3 = J(:,4);  
Jf4 = J(:,5); 
Jf5 = J(:,6);  
Jf6 = J(:,7); 
Jf7 = J(:,8);  
Jf8 = J(:,9); 
Jf9 = J(:,10); 
Jf10 = J(:,11);  
Jf11 = J(:,12); 
% Jf12 = J(:,13); 
% Jf13 = J(:,14); 

K = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_NNOutputs\Kdsyngain_ABOL_3(4).dat');   % Kd syn gain only
% K = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\KdsyngainABCLTest1(0.4).dat');   % Kd syn gain only
Kf1 = K(:,2);  
Kf2 = K(:,3);
Kf3 = K(:,4);  
Kf4 = K(:,5); 
Kf5 = K(:,6);  
Kf6 = K(:,7); 
Kf7 = K(:,8);  
Kf8 = K(:,9); 
Kf9 = K(:,10); 
Kf10 = K(:,11);  
Kf11 = K(:,12); 
% Kf12 = K(:,13); 
% Kf13 = K(:,14); 

L = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_NNOutputs\kdpregain_ABOL_3(4).dat');   % Kd syn gain only
% L = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\KdsyngainABCLTest1(0.4).dat');   % Kd syn gain only
Lf1 = L(:,2);  
Lf2 = L(:,3);
Lf3 = L(:,4);  
Lf4 = L(:,5); 
Lf5 = L(:,6);  
Lf6 = L(:,7); 
Lf7 = L(:,8);  
Lf8 = L(:,9); 
Lf9 = L(:,10); 
Lf10 = L(:,11);  
Lf11 = L(:,12); 
% Lf12 = L(:,13); 
% Lf13 = L(:,14); 

M = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_NNOutputs\ErrorScope_ABOL_3(4).dat');   % Kd syn gain only
% L = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\KdsyngainABCLTest1(0.4).dat');   % Kd syn gain only
Mf1 = M(:,2);  
Mf2 = M(:,3);
Mf3 = M(:,4);  
Mf4 = M(:,5); 
Mf5 = M(:,6);  
Mf6 = M(:,7); 
Mf7 = M(:,8);  
Mf8 = M(:,9); 
Mf9 = M(:,10); 
Mf10 = M(:,11);  
Mf11 = M(:,12); 
% Mf12 = M(:,13); 
% Mf13 = M(:,14); 

%% Plots

% Desired and Actual Position Figure
figure
plot(time,Bf6)  % system input
hold on
plot(time,Ef6)  % output
xlabel('Time')
ylabel('Amplitude (angle in nA)')
title('Frequency 1')
% legend('Desired Position','Actual Position','Location','northeast')
legend('Error','NNoutput','Location','northeast')
xlim([0 10])
ylim([-3 3])

%% System Input System Output for Closed Loop Figures
figure
plot(time,Af1)  % system input
hold on
% plot(time,Bf1)    % error
% plot(time,Ef1)      % NN output
plot(time,Ff1)      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 1')
% legend('Input Signal','Error','Location','northeast')
% legend('Error','NNoutput','Location','northeast')
legend('system Input','system Output','Location','northeast')
xlim([0 200])
ylim([-3 2.5])

% Desired and Actual Position Figure
figure
plot(time,Af2)  % system input
hold on
% plot(time,Bf2)    % error
% plot(time,Ef2)      % NN output
plot(time,Ff2)      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 2')
% legend('Input Signal','Error','Location','northeast')
% legend('Error','NNoutput','Location','northeast')
legend('system Input','system Output','Location','northeast')
xlim([0 100])
ylim([-3 2.5])

% Desired and Actual Position Figure
figure
plot(time,Af3)  % system input
hold on
plot(time,Bf3)    % error
% plot(time,Ef3)      % NN output
% plot(time,Ff3)      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 3')
% legend('Input Signal','Error','Location','northeast')
% legend('Error','NNoutput','Location','northeast')
legend('system Input','system Output','Location','northeast')
xlim([0 70])
ylim([-3 3])

% Desired and Actual Position Figure
figure
plot(time,Af4)  % system input
hold on
% plot(time,Bf4)    % error
% plot(time,Ef4)      % NN output
plot(time,Ff4)      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 4')
% legend('Input Signal','Error','Location','northeast')
% legend('Error','NNoutput','Location','northeast')
legend('system Input','system Output','Location','northeast')
xlim([0 40])
ylim([-20 20])

% Desired and Actual Position Figure
figure
plot(time,Af5)  % system input
hold on
% plot(time,Bf5)    % error
% plot(time,Ef5)      % NN output
plot(time,Ff5)      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 5')
% legend('Input Signal','Error','Location','northeast')
% legend('Error','NNoutput','Location','northeast')
legend('system Input','system Output','Location','northeast')
xlim([0 40])
ylim([-20 20])

% Desired and Actual Position Figure
figure
plot(time,Af6)  % system input
hold on
% plot(time,Bf6)    % error
% plot(time,Ef6)      % NN output
plot(time,Ff6)      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 6')
% legend('Input Signal','Error','Location','northeast')
% legend('Error','NNoutput','Location','northeast')
legend('system Input','system Output','Location','northeast')
xlim([0 30])
ylim([-20 20])


% Desired and Actual Position Figure
figure
plot(time,Af7)  % system input
hold on
% plot(time,Bf7)    % error
% plot(time,Ef7)      % NN output
plot(time,Ff7)      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 7')
% legend('Input Signal','Error','Location','northeast')
% legend('Error','NNoutput','Location','northeast')
legend('system Input','system Output','Location','northeast')
xlim([0 15])
ylim([-15 15])

% Desired and Actual Position Figure
figure
plot(time,Af8)  % system input
hold on
% plot(time,Bf8)    % error
% plot(time,Ef8)      % NN output
plot(time,Ff8)      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 8')
% legend('Input Signal','Error','Location','northeast')
% legend('Error','NNoutput','Location','northeast')
legend('system Input','system Output','Location','northeast')
xlim([0 10])
ylim([-10 10])

% Desired and Actual Position Figure
figure
plot(time,Af9)  % system input
hold on
% plot(time,Bf9)    % error
% plot(time,Ef9)      % NN output
plot(time,Ff9)      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 9')
% legend('Input Signal','Error','Location','northeast')
% legend('Error','NNoutput','Location','northeast')
legend('system Input','system Output','Location','northeast')
xlim([0 7])
ylim([-10 10])

% Desired and Actual Position Figure
figure
plot(time,Af10)  % system input
hold on
% plot(time,Bf10)    % error
% plot(time,Ef10)      % NN output
plot(time,Ff10)      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 10')
% legend('Input Signal','Error','Location','northeast')
% legend('Error','NNoutput','Location','northeast')
legend('system Input','system Output','Location','northeast')
xlim([0 5])
ylim([-5 5])

% Desired and Actual Position Figure
figure
plot(time,Af11)  % system input
hold on
% plot(time,Bf11)    % error
% plot(time,Ef11)      % NN output
plot(time,Ff11)      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 11')
% legend('Input Signal','Error','Location','northeast')
% legend('Error','NNoutput','Location','northeast')
legend('system Input','system Output','Location','northeast')
xlim([0 3])
ylim([-2 2])

% Desired and Actual Position Figure
figure
plot(time,Af12)  % system input
hold on
% plot(time,Bf12)    % error
% plot(time,Ef12)      % NN output
plot(time,Ff12)      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 12')
% legend('Input Signal','Error','Location','northeast')
% legend('Error','NNoutput','Location','northeast')
legend('system Input','system Output','Location','northeast')
xlim([0 2])
ylim([-1 1])
%% divide
figure
plot(time,Af12)  % system input
hold on
% plot(time,Bf12)    % error
% plot(time,Ef12)      % NN output
plot(time,Ff12)      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 13')
% legend('Input Signal','Error','Location','northeast')
% legend('Error','NNoutput','Location','northeast')
legend('system Input','system Output','Location','northeast')
xlim([0 2])
ylim([-1 1])

%% Error Input NN Output for Open Loop PD Figures
figure
hold on
plot(time,Bf1)   	% error
plot(time,Ef1)      % NN output
plot(time,Mf1)      % error scope
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 1')
legend('Error','NNoutput','Error Scope','Location','northeast')
xlim([0 70])
ylim([-5 5])

figure
hold on
plot(time,Bf2)   	% error
plot(time,Ef2)      % NN output
plot(time,Mf2)      % error scope
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 2')
legend('Error','NNoutput','Error Scope','Location','northeast')
xlim([0 60])
ylim([-5 5])

figure
hold on
plot(time,Bf3)   	% error
plot(time,Ef3)      % NN output
plot(time,Mf3)      % error scope
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 3')
legend('Error','NNoutput','Error Scope','Location','northeast')
xlim([0 60])
ylim([-5 5])

figure
hold on
plot(time,Bf4)   	% error
plot(time,Ef4)      % NN output
plot(time,Mf4)      % error scope
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 4')
legend('Error','NNoutput','Error Scope','Location','northeast')
xlim([0 60])
ylim([-5 5])

figure
hold on
plot(time,Bf5)   	% error
plot(time,Ef5)      % NN output
plot(time,Mf5)      % error scope
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 5')
legend('Error','NNoutput','Error Scope','Location','northeast')
xlim([0 50])
ylim([-5 5])

figure
hold on
plot(time,Bf6)   	% error
plot(time,Ef6)      % NN output
plot(time,Mf6)      % error scope
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 6')
legend('Error','NNoutput','Error Scope','Location','northeast')
xlim([0 12])
ylim([-20 20])


figure
hold on
plot(time,Bf7)   	% error
plot(time,Ef7)      % NN output
plot(time,Mf7)      % error scope
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 7')
legend('Error','NNoutput','Error Scope','Location','northeast')
xlim([0 10])
ylim([-20 20])


figure
hold on
plot(time,Bf8)   	% error
plot(time,Ef8)      % NN output
plot(time,Mf8)      % error scope
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 8')
legend('Error','NNoutput','Error Scope','Location','northeast')
xlim([0 8])
ylim([-20 20])

figure
hold on
plot(time,Bf9)   	% error
plot(time,Ef9)      % NN output
plot(time,Mf9)      % error scope
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 9')
legend('Error','NNoutput','Error Scope','Location','northeast')
xlim([0 7])
ylim([-20 20])

figure
hold on
plot(time,Bf10)   	% error
plot(time,Ef10)      % NN output
plot(time,Mf10)      % error scope
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 10')
legend('Error','NNoutput','Error Scope','Location','northeast')
xlim([0 5])
ylim([-20 20])

figure
hold on
plot(time,Bf11)   	% error
plot(time,Ef11)      % NN output
plot(time,Mf11)      % error scope
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 11')
legend('Error','NNoutput','Error Scope','Location','northeast')
xlim([0 3])
ylim([-20 20])

figure
hold on
plot(time,Bf12)   	% error
plot(time,Ef12)      % NN output
plot(time,Mf12)      % error scope
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 12')
legend('Error','NNoutput','Error Scope','Location','northeast')
xlim([0 2])
ylim([-1.5 1.5])
%% divide
figure
hold on
plot(time,Bf13)   	% error
plot(time,Ef13)      % NN output
plot(time,Mf13)      % error scope
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 13')
legend('Error','NNoutput','Error Scope','Location','northeast')
xlim([0 1])
ylim([-1.25 1.25])

%% Kp/Kd Figs

figure
hold on
plot(time,Bf1)      % error
plot(time,Cf1)      % Kp*Error
plot(time,Df1)      % Kd*de/dt
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 1')
legend('Error','Kp*Error','Kd*de/dt','Location','northeast')
% xlim([0 20]) % open loop
xlim([0 10])
ylim([-5 5])

figure
hold on
plot(time,Bf2)      % error
plot(time,Cf2)      % Kp*Error
plot(time,Df2)      % Kd*de/dt
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 2')
legend('Error','Kp*Error','Kd*de/dt','Location','northeast')
% xlim([0 10]) % open loop
xlim([0 10])
ylim([-5 5])

figure
hold on
plot(time,Bf3)      % error
plot(time,Cf3)      % Kp*Error
plot(time,Df3)      % Kd*de/dt
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 3')
legend('Error','Kp*Error','Kd*de/dt','Location','northeast')
% xlim([0 8]) % open loop
xlim([0 10])
ylim([-5 5])

figure
hold on
plot(time,Bf4)      % error
plot(time,Cf4)      % Kp*Error
plot(time,Df4)      % Kd*de/dt
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 4')
legend('Error','Kp*Error','Kd*de/dt','Location','northeast')
% xlim([0 6]) % open loop
xlim([0 10])
ylim([-5 5])
%% divide
figure
hold on
plot(time,Bf5)      % error
plot(time,Cf5)      % Kp*Error
plot(time,Df5)      % Kd*de/dt
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 5')
legend('Error','Kp*Error','Kd*de/dt','Location','northeast')
legend('Error','Kp*Error','Kd*de/dt','Location','northeast')
% xlim([0 4]) % open loop
xlim([0 30])
ylim([-3 4])

figure
hold on
plot(time,Bf6)      % error
plot(time,Cf6)      % Kp*Error
plot(time,Df6)      % Kd*de/dt
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 6')
legend('Error','Kp*Error','Kd*de/dt','Location','northeast')
% xlim([0 2]) % open loop
xlim([0 10])
ylim([-3 6])


figure
hold on
plot(time,Bf7)      % error
plot(time,Cf7)      % Kp*Error
plot(time,Df7)      % Kd*de/dt
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 7')
legend('Error','Kp*Error','Kd*de/dt','Location','northeast')
% xlim([0 0.6]) % open loop
xlim([0 10])
ylim([-2.25 5])

figure
hold on
plot(time,Bf8)      % error
plot(time,Cf8)      % Kp*Error
plot(time,Df8)      % Kd*de/dt
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 8')
legend('Error','Kp*Error','Kd*de/dt','Location','northeast')
% xlim([0 0.3]) % open loop
xlim([0 15])
ylim([-10 10])

figure
hold on
plot(time,Bf9)      % error
plot(time,Cf9)      % Kp*Error
plot(time,Df9)      % Kd*de/dt
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 9')
legend('Error','Kp*Error','Kd*de/dt','Location','northeast')
% xlim([0 0.2]) % open loop
xlim([0 10])
ylim([-10 10])

figure
hold on
plot(time,Bf10)      % error
plot(time,Cf10)      % Kp*Error
plot(time,Df10)      % Kd*de/dt
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 10')
legend('Error','Kp*Error','Kd*de/dt','Location','northeast')
xlim([0 2]) % open loop
% xlim([0 2])
ylim([-5 10])

figure
hold on
plot(time,Bf11)      % error
plot(time,Cf11)      % Kp*Error
plot(time,Df11)      % Kd*de/dt
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 11')
legend('Error','Kp*Error','Kd*de/dt','Location','northeast')
xlim([0 2])
ylim([-1 2])

figure
hold on
plot(time,Bf12)      % error
plot(time,Cf12)      % Kp*Error
plot(time,Df12)      % Kd*de/dt
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 12')
legend('Error','Kp*Error','Kd*de/dt','Location','northeast')
xlim([0 10])
ylim([-0.8 0.8])

%% System Input and Error and System Output
figure
plot(time,Af1)  % system input
hold on
plot(time,Bf1)    % error
plot(time,Ff1,':k')      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 1')
legend('Input Signal','Error','system output','Location','northeast')
xlim([0 200])
ylim([-4 4])

% Desired and Actual Position Figure
figure
plot(time,Af2)  % system input
hold on
plot(time,Bf2)    % error
plot(time,Ff2,':k')      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 2')
legend('Input Signal','Error','system output','Location','northeast')
xlim([0 180])
ylim([-4 4])

% Desired and Actual Position Figure
figure
plot(time,Af3)  % system input
hold on
plot(time,Bf3)    % error
plot(time,Ff3,':k')      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 3')
legend('Input Signal','Error','system output','Location','northeast')
xlim([0 120])
ylim([-4 4])

% Desired and Actual Position Figure
figure
plot(time,Af4)  % system input
hold on
plot(time,Bf4)    % error
plot(time,Ff4,':k')      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 4')
legend('Input Signal','Error','system output','Location','northeast')
xlim([0 80])
ylim([-2.5 2.5])

% Desired and Actual Position Figure
figure
plot(time,Af5)  % system input
hold on
plot(time,Bf5)    % error
plot(time,Ff5,':k')      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 5')
legend('Input Signal','Error','system output','Location','northeast')
xlim([0 40])
ylim([-3 3])

% Desired and Actual Position Figure
figure
plot(time,Af6)  % system input
hold on
plot(time,Bf6)    % error
plot(time,Ff6,':k')      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 6')
legend('Input Signal','Error','system output','Location','northeast')
xlim([0 20])
ylim([-3 2.5])

% Desired and Actual Position Figure
figure
plot(time,Af7)  % system input
hold on
plot(time,Bf7)    % error
plot(time,Ff7,':k')      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 7')
legend('Input Signal','Error','system output','Location','northeast')
xlim([0 10])
ylim([-2 2])

% Desired and Actual Position Figure
figure
plot(time,Af8)  % system input
hold on
plot(time,Bf8)    % error
plot(time,Ff8,':k')      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 8')
legend('Input Signal','Error','system output','Location','northeast')
xlim([0 10])
ylim([-2 2])

% Desired and Actual Position Figure
figure
plot(time,Af9)  % system input
hold on
plot(time,Bf9)    % error
plot(time,Ff9,':k')      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 9')
legend('Input Signal','Error','system output','Location','northeast')
xlim([0 100])
ylim([-3 2.5])

% Desired and Actual Position Figure
figure
plot(time,Af10)  % system input
hold on
plot(time,Bf10)    % error
plot(time,Ff10,':k')      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 10')
legend('Input Signal','Error','system output','Location','northeast')
xlim([0 100])
ylim([-3 2.5])

% Desired and Actual Position Figure
figure
plot(time,Af11)  % system input
hold on
plot(time,Bf11)    % error
plot(time,Ff11,':k')      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 11')
legend('Input Signal','Error','system output','Location','northeast')
xlim([0 100])
ylim([-3 2.5])
%% divide
% Desired and Actual Position Figure
figure
plot(time,Af12)  % system input
hold on
plot(time,Bf12)    % error
plot(time,Ff12,':k')      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 12')
legend('Input Signal','Error','system output','Location','northeast')
xlim([0 100])
ylim([-3 2.5])

% Desired and Actual Position Figure
figure
plot(time,Af13)  % system input
hold on
plot(time,Bf13)    % error
plot(time,Ff13,':k')      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Frequency 13')
legend('Input Signal','Error','system output','Location','northeast')
xlim([0 100])
ylim([-3 2.5])