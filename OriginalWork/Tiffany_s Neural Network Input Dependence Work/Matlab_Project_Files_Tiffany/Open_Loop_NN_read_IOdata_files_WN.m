% Read data from NNinput and output data files

% MAKE SURE TO UPDATE FILENAMES AND AMPLITUDES IN ()

clear
clc

% Import Data
% assumes all have the same time scale
A = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\NNinput_WNOL_61(0.2).dat');   % Desired input to NN
% A = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\NNinputABCLTest1(0.4).dat');   % Desired input to NN
time = A(:,1);      % s
Af1 = A(:,2);  

B = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\ErrorN3_WNOL_61(0.2).dat');   % Error N3 outputs
% B = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\ErrorN3ABCLTest1(0.4).dat');   % Error N3 outputs
Bf1 = B(:,2);  

C = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\Kp_WNOL_61(0.2).dat');   % Kp*Error circuit outputs
% C = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\KpOutputABCLTest1(0.4).dat');   % Kp*Error circuit outputs
Cf1 = C(:,2);  

D = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\Kd_WNOL_61(0.2).dat');   % Kd*de/dt circuit outputs
% D = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\KdOutputABCLTest1(0.4).dat');   % Kd*de/dt circuit outputs
Df1 = D(:,2);  

E = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\NNOutput_WNOL_61(0.2).dat');   % NN circuit outputs
% E = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\NNOutputABCLTest1(0.4).dat');   % NN circuit outputs
Ef1 = E(:,2);  

% F = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\NNangleOut_WNCL_34(0.6).dat');   % Actual output (input to NN)
% % F = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\NNangleOutABCLTest1(0.4).dat');   % Actual output (input to NN)
% Ff1 = F(:,2);  

G = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\Kpscaler_WNOL_61(0.2).dat');   % Kp scaler
% G = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\KpscalerABCLTest1(0.4).dat');   % Kp scaler
Gf1 = G(:,2);  

H = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\Kdscaler_WNOL_61(0.2).dat');   % Kd scaler
% H = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\KdscalerABCLTest1(0.4).dat');   % Kd scaler
Hf1 = H(:,2);  

J = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\Kpsyngain_WNOL_61(0.2).dat');   % Kp syn gain only
% J = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\KpsyngainABCLTest1(0.4).dat');   % Kp syn gain only
Jf1 = J(:,2);  

K = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\Kdsyngain_WNOL_61(0.2).dat');   % Kd syn gain only
% K = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\KdsyngainABCLTest1(0.4).dat');   % Kd syn gain only
Kf1 = K(:,2);  

L = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\kdpregain_WNOL_61(0.2).dat');   % Kd syn gain only
% L = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\KdsyngainABCLTest1(0.4).dat');   % Kd syn gain only
Lf1 = L(:,2);  

M = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\ErrorScope_WNOL_61(0.2).dat');   % Kd syn gain only
% L = importdata('C:\Users\Tiffany\Desktop\B_NN_Outputs\KdsyngainABCLTest1(0.4).dat');   % Kd syn gain only
Mf1 = M(:,2);  

%% Plots

% Desired and Actual Position Figure
figure
plot(time,Af1)  % system input
hold on
plot(time,Ef1)  % output
xlabel('Time (s)')
ylabel('Amplitude (nA)')
title('Whitenoise Input/Output')
% legend('Desired Position','Actual Position','Location','northeast')
legend('Error','NNoutput','Location','northeast')
xlim([0 5])
% ylim([-20 20])
set(gcf, 'Color', 'w');


% System Input System Output for Closed Loop Figures
figure
plot(time,Af1)  % system input
hold on
% plot(time,Bf1)    % error
% plot(time,Ef1)      % NN output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Whitenoise Test')
% legend('Input Signal','Error','Location','northeast')
% legend('Error','NNoutput','Location','northeast')
legend('system Input','Location','northeast')
% xlim([0 200])
% ylim([-3 2.5])

mWindow = 1000;
errorMean = movmean(abs(Bf1),mWindow);
torqueMean = movmean(abs(Ef1),mWindow);
%%
% Error Input NN Output for Open Loop PD Figures
figure
hold on
plot(time,Bf1)   	% error
plot(time,Ef1)      % NN output
% plot(time,Mf1)      % error scope
plot(time,errorMean,'Linewidth',1.4)
plot(time,errorTorque,'Linewidth',1.4)
xlabel('Time')
ylabel('Amplitude (nA)')
title('Whitenoise Test')
% legend('Error','NNoutput','Error Scope','Location','northeast')
legend('Error','NNoutput','Error Mean','Torque Mean','Location','northeast')
% xlim([0 70])
% ylim([-5 5])
%%
% Kp/Kd Figs

figure
hold on
plot(time,Bf1)      % error
plot(time,Cf1)      % Kp*Error
plot(time,Df1)      % Kd*de/dt
xlabel('Time')
ylabel('Amplitude (nA)')
title('Whitenoise Test')
legend('Error','Kp*Error','Kd*de/dt','Location','northeast')
% xlim([0 20]) % open loop
% xlim([0 75])
% ylim([-3 3])

% System Input and Error and System Output
figure
plot(time,Af1)  % system input
hold on
plot(time,Bf1)    % error
% plot(time,Ff1,':k')      % system output
xlabel('Time')
ylabel('Amplitude (nA)')
title('Whitenoise Test')
legend('Input Signal','Error','Location','northeast')
% xlim([0 200])
% ylim([-4 4])
