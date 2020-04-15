% Compare bodes for closed loop Controller overall, Kp and Kd circuits for varying
% inputs

clear
clc

% Import Data
% assumes all have the same frequency scale
A = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\SYS_ABCL_9(0.1).dat');   % sys circuit outputs
freq = A(:,1);      % rad/s
magRatioCL_a = A(:,2);  
phaseCL_a = A(:,3);
CLstd_a = A(:,4);

P2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\PD_ABCL_9(0.1).dat');   % PD circuit
magRatioPD_a = P2(:,2);  
phasePD_a = P2(:,3);

B = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kp_ABCL_9(0.1).dat');   % Kp circuit outputs
magRatioKp_a = B(:,2);  
phaseKp_a = B(:,3);
error_a = B(:,4);

C = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kd_ABCL_9(0.1).dat');   % Kd circuit outputs
magRatioKd_a = C(:,2);
phaseKd_a = C(:,3);

F1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kdsyngain_ABCL_9(0.1).dat');   % Kd circuit outputs
Kdsyngain_a = F1(:,2);
Kdstd_a = F1(:,4);

Z2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\kdbase_ABCL_9(0.1).dat');   % Kp syn gain
kdbasegain_a = Z2(:,2);

P1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kpsyngain_ABCL_9(0.1).dat');   % Kp syn gain
Kpsyngain_a = P1(:,2);
Kpstd_a = P1(:,4);

D = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\SYS_ABCL_9(0.25).dat');   % sys circuit outputs
magRatioCL_b = D(:,2);  
phaseCL_b = D(:,3);

Q2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\PD_ABCL_9(0.25).dat');   % PD circuit
magRatioPD_b = Q2(:,2);  
phasePD_b = Q2(:,3);

E = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kp_ABCL_9(0.25).dat');   % Kp circuit outputs
magRatioKp_b = E(:,2);  
phaseKp_b = E(:,3);
error_b = E(:,4);

F = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kd_ABCL_9(0.25).dat');   % Kd circuit outputs
magRatioKd_b = F(:,2);
phaseKd_b = F(:,3);

A3 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\kdbase_ABCL_9(0.25).dat');   % Kp syn gain
kdbasegain_b = A3(:,2);

G1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kdsyngain_ABCL_9(0.25).dat');   % Kd circuit outputs
Kdsyngain_b = G1(:,2);

Q1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kpsyngain_ABCL_9(0.25).dat');   % Kp syn gain
Kpsyngain_b = Q1(:,2);

G = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\SYS_ABCL_9(0.5).dat');   % sys circuit outputs
magRatioCL_c = G(:,2);  
phaseCL_c = G(:,3);

R2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\PD_ABCL_9(0.5).dat');   % PD circuit
magRatioPD_c = R2(:,2);  
phasePD_c = R2(:,3);

H = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kp_ABCL_9(0.5).dat');   % Kp circuit outputs
magRatioKp_c = H(:,2);  
phaseKp_c = H(:,3);
error_c = H(:,4);

J = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kd_ABCL_9(0.5).dat');   % Kd circuit outputs
magRatioKd_c = J(:,2);
phaseKd_c = J(:,3);

B3 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\kdbase_ABCL_9(0.5).dat');   % Kp syn gain
kdbasegain_c = B3(:,2);

H1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kdsyngain_ABCL_9(0.5).dat');   % Kd circuit outputs
Kdsyngain_c = H1(:,2);

R1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kpsyngain_ABCL_9(0.5).dat');   % Kp syn gain
Kpsyngain_c = R1(:,2);

K = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\SYS_ABCL_9(0.75).dat');   % sys circuit outputs
magRatioCL_d = K(:,2);  
phaseCL_d = K(:,3);

T2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\PD_ABCL_9(0.75).dat');   % PD circuit
magRatioPD_d = T2(:,2);  
phasePD_d = T2(:,3);

L = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kp_ABCL_9(0.75).dat');   % Kp circuit outputs
magRatioKp_d = L(:,2);  
phaseKp_d = L(:,3);
error_d = L(:,4);

P = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kd_ABCL_9(0.75).dat');   % Kd circuit outputs
magRatioKd_d = P(:,2);
phaseKd_d = P(:,3);

C3 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\kdbase_ABCL_9(0.75).dat');   % Kp syn gain
kdbasegain_d = C3(:,2);

J1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kdsyngain_ABCL_9(0.75).dat');   % Kd circuit outputs
Kdsyngain_d = J1(:,2);

S1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kpsyngain_ABCL_9(0.75).dat');   % Kp syn gain
Kpsyngain_d = S1(:,2);

Q = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\SYS_ABCL_9(1).dat');   % sys circuit outputs
magRatioCL_e = Q(:,2);  
phaseCL_e = Q(:,3);

U2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\PD_ABCL_9(1).dat');   % PD circuit
magRatioPD_e = U2(:,2);  
phasePD_e = U2(:,3);

R = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kp_ABCL_9(1).dat');   % Kp circuit outputs
magRatioKp_e = R(:,2);  
phaseKp_e = R(:,3);
error_e = R(:,4);

T = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kd_ABCL_9(1).dat');   % Kd circuit outputs
magRatioKd_e = T(:,2);
phaseKd_e = T(:,3);

D3 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\kdbase_ABCL_9(1).dat');   % Kp syn gain
kdbasegain_e = D3(:,2);

K1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kdsyngain_ABCL_9(1).dat');   % Kd circuit outputs
Kdsyngain_e = K1(:,2);

T1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kpsyngain_ABCL_9(1).dat');   % Kp syn gain
Kpsyngain_e = T1(:,2);

U = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\SYS_ABCL_9(1.5).dat');   % sys circuit outputs
magRatioCL_f = U(:,2);  
phaseCL_f = U(:,3);

V2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\PD_ABCL_9(1.5).dat');   % PD circuit
magRatioPD_f = V2(:,2);  
phasePD_f = V2(:,3);

V = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kp_ABCL_9(1.5).dat');   % Kp circuit outputs
magRatioKp_f = V(:,2);  
phaseKp_f = V(:,3);
error_f = V(:,4);

W = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kd_ABCL_9(1.5).dat');   % Kd circuit outputs
magRatioKd_f = W(:,2);
phaseKd_f = W(:,3);

E3 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\kdbase_ABCL_9(1.5).dat');   % Kp syn gain
kdbasegain_f = E3(:,2);

L1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kdsyngain_ABCL_9(1.5).dat');   % Kd circuit outputs
Kdsyngain_f = L1(:,2);

U1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kpsyngain_ABCL_9(1.5).dat');   % Kp syn gain
Kpsyngain_f = U1(:,2);

X = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\SYS_ABCL_9(2).dat');   % sys circuit outputs
magRatioCL_g = X(:,2);  
phaseCL_g = X(:,3);

W2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\PD_ABCL_9(2).dat');   % PD circuit
magRatioPD_g = W2(:,2);  
phasePD_g = W2(:,3);

Y = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kp_ABCL_9(2).dat');   % Kp circuit outputs
magRatioKp_g = Y(:,2);  
phaseKp_g = Y(:,3);
error_g = Y(:,4);

Z = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kd_ABCL_9(2).dat');   % Kd circuit outputs
magRatioKd_g = Z(:,2);
phaseKd_g = Z(:,3);

F3 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\kdbase_ABCL_9(2).dat');   % Kp syn gain
kdbasegain_g = F3(:,2);

M1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kdsyngain_ABCL_9(2).dat');   % Kd circuit outputs
Kdsyngain_g = M1(:,2);

V1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kpsyngain_ABCL_9(2).dat');   % Kp syn gain
Kpsyngain_g = V1(:,2);

A1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\SYS_ABCL_9(3).dat');   % sys circuit outputs
magRatioCL_h = A1(:,2);  
phaseCL_h = A1(:,3);

X2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\PD_ABCL_9(3).dat');   % PD circuit
magRatioPD_h = X2(:,2);  
phasePD_h = X2(:,3);

B1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kp_ABCL_9(3).dat');   % Kp circuit outputs
magRatioKp_h = B1(:,2);  
phaseKp_h = B1(:,3);
error_h = B1(:,4);

C1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kd_ABCL_9(3).dat');   % Kd circuit outputs
magRatioKd_h = C1(:,2);
phaseKd_h = C1(:,3);

N1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kdsyngain_ABCL_9(3).dat');   % Kd circuit outputs
Kdsyngain_h = N1(:,2);

G3 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\kdbase_ABCL_9(3).dat');   % Kp syn gain
kdbasegain_h = G3(:,2);

W1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kpsyngain_ABCL_9(3).dat');   % Kp syn gain
Kpsyngain_h = W1(:,2);

J2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\SYS_ABCL_9(4).dat');   % sys circuit outputs
magRatioCL_j = J2(:,2);  
phaseCL_j = J2(:,3);

Y2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\PD_ABCL_9(4).dat');   % PD circuit
magRatioPD_j = Y2(:,2);  
phasePD_j = Y2(:,3);

K2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kp_ABCL_9(4).dat');   % Kp circuit outputs
magRatioKp_j = K2(:,2);  
phaseKp_j = K2(:,3);
error_j = K2(:,4);
CLstd_j = A(:,4);

L2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kd_ABCL_9(4).dat');   % Kd circuit outputs
magRatioKd_j = L2(:,2);
phaseKd_j = L2(:,3);

H3 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\kdbase_ABCL_9(4).dat');   % Kp syn gain
kdbasegain_j = H3(:,2);

M2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kdsyngain_ABCL_9(4).dat');   % Kd circuit outputs
Kdsyngain_j = M2(:,2);
Kdstd_j = M2(:,4);

N2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABCL_Outputs\Kpsyngain_ABCL_9(4).dat');   % Kp syn gain
Kpsyngain_j = N2(:,2);
Kpstd_j = N2(:,4);

% Closed Loop White Noise Data
D1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\NNinput_WNCL_32(0.6).dat');   % Desired input to NN
timeWNCL = D1(:,1);      % s
NNinput_desWNCL = D1(:,2);  

E1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\NNangleOut_WNCL_32(0.6).dat');   
SysOutputWNCL = E1(:,2);  

A2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\ErrorScope_WNCL_32(0.6).dat');   
ErrorScopeWNCL = A2(:,2);  

B2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\Kp_WNCL_32(0.6).dat');   
KpScopeWNCL = B2(:,2); 

C2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\Kd_WNCL_32(0.6).dat');   
KdScopeWNCL = C2(:,2); 

D2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\kdpregain_WNCL_32(0.6).dat');   
kdpregainWNCL = D2(:,2); 

E2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\Kpscaler_WNCL_32(0.6).dat');   
KpscalerWNCL = E2(:,2); 

F2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\Kdscaler_WNCL_32(0.6).dat');   
KdscalerWNCL = F2(:,2); 

G2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\NNOutput_WNCL_32(0.6).dat');   
NNOutputWNCL = G2(:,2); 

H2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\ErrorN3_WNCL_32(0.6).dat');   
ErrorN3WNCL = H2(:,2); 

freq_Hz = freq/(2*pi);
freq_Hz_flip = [freq_Hz,fliplr(freq_Hz)];

%% Process White Noise Data

dtWNCL = timeWNCL(3,1)-timeWNCL(2,1);
TsWNCL = 1/dtWNCL;
numSegmentsWNCL = 4;   % # segments to divide signal into for fft averaging
movAvgWinTimeWNCL = 0.5;   % window size for moving average in seconds

[f_fftWNCL,magSys,phaseSys,AmpinSys,AmpoutSys] = bodebyFFT(NNinput_desWNCL,SysOutputWNCL,TsWNCL,numSegmentsWNCL,movAvgWinTimeWNCL);
[f_fftWNCL,magKpgain,phaseKpgain,AmpinKpgain,AmpoutKpgain] = bodebyFFT(ErrorScopeWNCL,KpScopeWNCL,TsWNCL,numSegmentsWNCL,movAvgWinTimeWNCL);
[f_fftWNCL,magKdgain,phaseKdgain,AmpinKdgain,AmpoutKdgain] = bodebyFFT(kdpregainWNCL,KdScopeWNCL,TsWNCL,numSegmentsWNCL,movAvgWinTimeWNCL);
[f_fftWNCL,magkdbase,phasekdbase,Ampinkdbase,Ampoutkdbase] = bodebyFFT(ErrorScopeWNCL,kdpregainWNCL,TsWNCL,numSegmentsWNCL,movAvgWinTimeWNCL);
[f_fftWNCL,magPD,phasePD,AmpinPD,AmpoutPD] = bodebyFFT(ErrorN3WNCL,NNOutputWNCL,TsWNCL,numSegmentsWNCL,movAvgWinTimeWNCL);

gainHuman = 833;
gainAB = 833;
gainWN = 1388;
tau = 0.1;

Kpstart = 0.1; Kpend = 2.5; Kpdivs = 40;
Kdstart = 0.1; Kdend = 1; Kddivs = 40;
[KpFitWN,KdFitWN] = LSRfitCLPD_MagPhase_2Param(Kpstart,Kpend,Kpdivs,Kdstart,Kdend,Kddivs,f_fftWNCL(1:floor(end/50)),tau,gainWN,magSys(1:floor(end/50)),phaseSys(1:floor(end/50)))
%%
mWindow_freq = 50;
mWindow_time = 5000;
Kpmean = movmean(magKpgain(1:floor(end/50)),mWindow_freq);
Kpavgmean = mean(Kpmean)
Kpstd = movstd(magKpgain(1:floor(end/50)),mWindow_freq);
Kpavgstd = mean(Kpstd)

kdbasemean = movmean(magkdbase(1:floor(end/50)),mWindow_freq);
kdbaseavgmean = mean(kdbasemean)
kdbasestd = movstd(magkdbase(1:floor(end/50)),mWindow_freq);
kdbaseavgstd = mean(kdbasestd)

Kdmean = movmean(magKdgain(1:floor(end/50)),mWindow_freq)*kdbaseavgmean;
Kdavgmean = mean(Kdmean)*kdbaseavgmean
Kdstd = movstd(magKdgain(1:floor(end/50)),mWindow_freq)*kdbaseavgmean;
Kdavgstd = mean(Kdstd)*kdbaseavgmean
KpscalerMean = movmean(KpscalerWNCL,mWindow_time);
KdscalerMean = movmean(KdscalerWNCL,mWindow_time);
ErrorMean = movmean(ErrorScopeWNCL,mWindow_time);
KpScopeMean = movmean(KpScopeWNCL,mWindow_time);
NNOutputMean = movmean(abs(NNOutputWNCL),mWindow_time);
ErrorN3Mean = movmean(abs(ErrorN3WNCL),mWindow_time);

figure
semilogx(f_fftWNCL(1:floor(end/50)),Kpmean,'Linewidth',1.4)
hold on
semilogx(f_fftWNCL(1:floor(end/50)),Kpmean+Kpstd,'--','Linewidth',1.4)
semilogx(f_fftWNCL(1:floor(end/50)),max(Kpmean-Kpstd,0),'--','Linewidth',1.4)
semilogx(f_fftWNCL(1:floor(end/50)),magKpgain(1:floor(end/50)),'o')
legend('Kp mean','Kp + 1stdev','Kp - 1stdev','data')
xlim([0.01 10])

figure
semilogx(f_fftWNCL(1:floor(end/50)),Kdmean,'Linewidth',1.4)
hold on
semilogx(f_fftWNCL(1:floor(end/50)),Kdmean+Kdstd,'--','Linewidth',1.4)
semilogx(f_fftWNCL(1:floor(end/50)),max(Kdmean-Kdstd,0),'--','Linewidth',1.4)
semilogx(f_fftWNCL(1:floor(end/50)),magKdgain(1:floor(end/50)),'o')
legend('Kd mean','Kd + 1stdev','Kd - 1stdev','data')
xlim([0.01 10])

figure
semilogx(f_fftWNCL(1:floor(end/50)),kdbasemean,'Linewidth',1.4)
hold on
semilogx(f_fftWNCL(1:floor(end/50)),kdbasemean + kdbasestd,'--','Linewidth',1.4)
semilogx(f_fftWNCL(1:floor(end/50)),max(kdbasemean - kdbasestd,0),'--','Linewidth',1.4)
semilogx(f_fftWNCL(1:floor(end/50)),magkdbase(1:floor(end/50)),'o')
legend('kd base mean','kd base + 1stdev','kd base - 1stdev','data')
xlim([0.01 10])

figure
plot(timeWNCL,KpscalerWNCL,'o')
hold on
plot(timeWNCL,KpscalerMean,'Linewidth',1.4)
title('Kp Scaler')

figure
plot(timeWNCL,KdscalerWNCL,'o')
hold on
plot(timeWNCL,KdscalerMean,'Linewidth',1.4)
title('Kd Scaler')

figure
plot(timeWNCL,ErrorScopeWNCL)
hold on
plot(timeWNCL,ErrorMean,'Linewidth',2)
plot(timeWNCL,KpScopeWNCL)
plot(timeWNCL,KpScopeMean,'Linewidth',2)
title('Kp Gain over Time')
legend('Error Scope','Error Mean','Kp Scope','Kp Mean')

figure
plot(timeWNCL,ErrorN3WNCL)
hold on
plot(timeWNCL,ErrorN3Mean,'Linewidth',2)
plot(timeWNCL,NNOutputWNCL)
plot(timeWNCL,NNOutputMean,'Linewidth',2)
title('Error and Torque over Time')
legend('Error','Error Mean','Output Torque','Torque Mean')

%% Closed Loop Human Reference Model & Auto Bode Fits
s = tf('s');

% Peterka's "normal subject" constants
m = 83.3; %kg 
h = .896; %m - center of mass height
I = 81.1; %kg m^2
g = 9.81;% m/s

Gp = 1/(I*s^2-m*g*h);

% Fit TF to Human Data (STATIC)
% gain = 833;
tauHuman = 0.1;
Kd = 0.05*6;   % for gain = 1, Kd = 250
Kp = 0.21*6;   % for gain = 1, Kp = 1050

% fb1 = 5;         % cutoff frequency in Hz (base 2)
% lowFiltDerivative = (1/(s/(2*pi*fb1)+1));  % low pass filter
% 
% fb2 = 80;         % cutoff frequency in Hz (base 80)
% lowFiltNeuron = (1/(s/(2*pi*fb2)+1));  % low pass filter

[numHuman,denHuman] = pade(tauHuman,1);    % pade delay approximation first order (in order to use impulse command)
delayHuman = tf(numHuman,denHuman);

[num,den] = pade(tau,1);    % pade delay approximation first order (in order to use impulse command)
delay = tf(num,den);

Gc = Kp + Kd*s;
oltf = delayHuman*Gc*Gp*gainHuman;

cltf = feedback(oltf,1);

[magHuman,phaseHuman,omegaHuman] = bode(cltf);    % requires frequency in rad/s
magHuman = squeeze(magHuman);
phaseHuman = squeeze(phaseHuman);
f_Human = omegaHuman/(2*pi);

% Closed loop fits to Auto Bode Data

Kpstart = 0.5; Kpend = 2.5; Kpdivs = 40;
Kdstart = 0.1; Kdend = 0.7; Kddivs = 40;
[KpFitLow,KdFitLow] = LSRfitCLPD_MagPhase_2Param(Kpstart,Kpend,Kpdivs,Kdstart,Kdend,Kddivs,freq_Hz,tau,gainAB,magRatioCL_a,phaseCL_a)
[KpFitHigh,KdFitHigh] = LSRfitCLPD_MagPhase_2Param(Kpstart,Kpend,Kpdivs,Kdstart,Kdend,Kddivs,freq_Hz,tau,gainAB,magRatioCL_j,phaseCL_j)


GcHigh = KpFitHigh + KdFitHigh*s;
oltfHigh = delay*GcHigh*Gp*gainAB;
cltfHigh = feedback(oltfHigh,1);

[magHigh,phaseHigh] = bode(cltfHigh,freq);    % requires frequency in rad/s
magHigh = squeeze(magHigh);
phaseHigh = squeeze(phaseHigh);

GcLow = KpFitLow + KdFitLow*s;
oltfLow = delay*GcLow*Gp*gainAB;
cltfLow = feedback(oltfLow,1);

[magLow,phaseLow] = bode(cltfLow,freq);    % requires frequency in rad/s
magLow = squeeze(magLow);
phaseLow = squeeze(phaseLow);
%%
% magflip = [magHigh,fliplr(magLow)];
% 
% patch(freq_Hz_flip,magflip,'b','FaceAlpha',0.1)
% figure
% hold on
% area(freq_Hz,magLow,'FaceColor','b')
% area(freq_Hz,magHigh,'FaceColor','w')
%%
% Closed loop fit to whitenoise data
% KpFitWN = 1;
% KdFitWN = 0.25;

% KpFitWN = 1;
% KdFitWN = 0.2;

GcWN = KpFitWN + KdFitWN*s;
oltfWN = delay*GcWN*Gp*gainWN;
cltfWN = feedback(oltfWN,1);

f_fftWNCL_rad = f_fftWNCL*2*pi;
[magWN,phaseWN] = bode(cltfWN,f_fftWNCL_rad);    % requires frequency in rad/s
magWN = squeeze(magWN);
phaseWN = squeeze(phaseWN);

%% Plots
% PD Open loop figure
figure
subplot(2,1,1)
semilogx(f_fftWNCL,magPD,'o')
hold on
% semilogx(f_fftWNCL,magWN,'Linewidth',1.2)
semilogx(freq_Hz,magRatioPD_a,'-s','Linewidth',1.4)
semilogx(freq_Hz,magRatioPD_b,'-+','Linewidth',1.4)
semilogx(freq_Hz,magRatioPD_c,'-^','Linewidth',1.4)
semilogx(freq_Hz,magRatioPD_d,'-p','Linewidth',1.4)
semilogx(freq_Hz,magRatioPD_e,'-*','Linewidth',1.4)
semilogx(freq_Hz,magRatioPD_f,'-x','Linewidth',1.4)
semilogx(freq_Hz,magRatioPD_g,'-s','Linewidth',1.4)
semilogx(freq_Hz,magRatioPD_h,'-h','Linewidth',1.4)
semilogx(freq_Hz,magRatioPD_j,'-d','Linewidth',1.4)
% semilogx(f_Human,magHuman,'Linewidth',1.2)
% semilogx(freq_Hz,magHigh,'--','Linewidth',1.4)
% semilogx(freq_Hz,magLow,'--','Linewidth',1.4)
xlim([0.1 2.25])
% fill(freq_Hz_flip,magflip,'b','FaceAlpha',0.1)
xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('Open Loop PD Transfer Function with Variable Inputs')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','3nA input','White Noise Data','Human Fit TF','Auto Bode Fit High','Auto Bode Fit Low','Location','southwest')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','White Noise Data','Human Fit TF','Auto Bode Fit High','Auto Bode Fit Low','Location','northeast')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','White Noise Data','Human Fit TF','Auto Bode Fit High','Auto Bode Fit Low','Location','northwest')
% legend('0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','White Noise Data','Human Fit TF','Auto Bode Fit High','Auto Bode Fit Low','Location','northeast')
% legend('White Noise Data','White Noise Fit','0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','3nA input','4nA input','Human Fit TF','Auto Bode Fit High','Auto Bode Fit Low','Location','northwest')
legend('White Noise Data','0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','3nA input','4nA input','Location','northwest')

% Error figure
subplot(2,1,2)
semilogx(freq_Hz,error_a,'-o')
hold on
semilogx(freq_Hz,error_b,'-+')

semilogx(freq_Hz,error_c,'-^')

semilogx(freq_Hz,error_d,'-p')
semilogx(freq_Hz,error_e,'-*')
semilogx(freq_Hz,error_f,'-x')
semilogx(freq_Hz,error_g,'-s')
semilogx(freq_Hz,error_h,'-h')
semilogx(freq_Hz,error_j,'-d')
xlabel('Frequency (Hz)')
ylabel('Average Error Signal Amplitude')
title('Average Error Amplitudes')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','3nA input','Location','northwest')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','Location','northwest')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','Location','northwest')
% legend('0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','Location','northwest')
legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','3nA input','4nA input','Location','northwest')
xlim([0.1 2.25])

% ylim([-40 20])
% xlim([0.01 10])
set(gcf, 'Color', 'w');

% CLosed Loop figure
figure
subplot(2,1,1)
% semilogx(f_fftWNCL,magSys,'o')
% semilogx(f_fftWNCL,magWN,'Linewidth',1.4)
semilogx(freq_Hz,magRatioCL_a,'-o','Linewidth',1.2)
hold on
semilogx(freq_Hz,magRatioCL_b,'-+','Linewidth',1.2)
semilogx(freq_Hz,magRatioCL_c,'-^','Linewidth',1.2)
semilogx(freq_Hz,magRatioCL_d,'-p','Linewidth',1.2)
semilogx(freq_Hz,magRatioCL_e,'-*','Linewidth',1.2)
semilogx(freq_Hz,magRatioCL_f,'-x','Linewidth',1.2)
semilogx(freq_Hz,magRatioCL_g,'-s','Linewidth',1.2)
semilogx(freq_Hz,magRatioCL_h,'-h','Linewidth',1.2)
semilogx(freq_Hz,magRatioCL_j,'-d','Linewidth',1.2)
% semilogx(f_Human,magHuman,'Linewidth',1.2)
% semilogx(freq_Hz,magHigh,'--','Linewidth',1.4)
% semilogx(freq_Hz,magLow,'--','Linewidth',1.4)
xlim([0.1 2.25])
% fill(freq_Hz_flip,magflip,'b','FaceAlpha',0.1)
xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('Closed Loop System Magnitude with Variable Inputs')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','3nA input','White Noise Data','Human Fit TF','Auto Bode Fit High','Auto Bode Fit Low','Location','southwest')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','White Noise Data','Human Fit TF','Auto Bode Fit High','Auto Bode Fit Low','Location','northeast')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','White Noise Data','Human Fit TF','Auto Bode Fit High','Auto Bode Fit Low','Location','northwest')
% legend('0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','White Noise Data','Human Fit TF','Auto Bode Fit High','Auto Bode Fit Low','Location','northeast')
% legend('White Noise Data','White Noise Fit','0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','3nA input','4nA input','Location','northwest')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','3nA input','4nA input','Location','northwest')


% ylim([-40 20])
% xlim([0.01 10])
set(gcf, 'Color', 'w');
% export_fig test.bmp -q101
% export_fig test.png -q101
% export_fig('C:\Users\Tiffany\Desktop\Sim_Outputs\Figs\CLABWN_MultiInput_Weighting1_WN26', '-bmp','-q101','-nocrop','-a1');

% Error figure
subplot(2,1,2)
semilogx(freq_Hz,error_a,'-o')
hold on
semilogx(freq_Hz,error_b,'-+')

semilogx(freq_Hz,error_c,'-^')

semilogx(freq_Hz,error_d,'-p')
semilogx(freq_Hz,error_e,'-*')
semilogx(freq_Hz,error_f,'-x')
semilogx(freq_Hz,error_g,'-s')
semilogx(freq_Hz,error_h,'-h')
semilogx(freq_Hz,error_j,'-d')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Average Error Signal Amplitudes')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','3nA input','Location','northwest')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','Location','northwest')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','Location','northwest')
% legend('0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','Location','northwest')
legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','3nA input','4nA input','Location','northwest')
xlim([0.1 2.25])

% Kd figure
figure
subplot(2,1,1)
semilogx(freq_Hz,magRatioKd_a,'-o')
hold on
semilogx(freq_Hz,magRatioKd_b,'-+')

semilogx(freq_Hz,magRatioKd_c,'-^')
semilogx(freq_Hz,magRatioKd_d,'-p')
semilogx(freq_Hz,magRatioKd_e,'-*')
semilogx(freq_Hz,magRatioKd_f,'-x')
semilogx(freq_Hz,magRatioKd_g,'-s')
semilogx(freq_Hz,magRatioKd_h,'-h')
semilogx(freq_Hz,magRatioKd_j,'-d')
% semilogx(freq_Hz,magKd)
xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('Kd Circuit Transfer Function Variable Inputs')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','3nA input','Previous fit for 1nA','Location','southwest')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','Location','southwest')
% ylim([-40 20])
set(gcf, 'Color', 'w');

% Kp figure
subplot(2,1,2)
semilogx(freq_Hz,magRatioKp_a,'-o')
hold on
semilogx(freq_Hz,magRatioKp_b,'-+')

semilogx(freq_Hz,magRatioKp_c,'-^')

semilogx(freq_Hz,magRatioKp_d,'-p')
semilogx(freq_Hz,magRatioKp_e,'-*')
semilogx(freq_Hz,magRatioKp_f,'-x')
semilogx(freq_Hz,magRatioKp_g,'-s')
semilogx(freq_Hz,magRatioKp_h,'-h')
semilogx(freq_Hz,magRatioKp_j,'-d')
% semilogx(freq_Hz,magKp)
xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('Kp Circuit Transfer Function Variable Inputs')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','3nA input','Location','southwest')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','Location','southwest')
% legend('0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','Location','southwest')
legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','3nA input','4nA input','Location','southwest')

% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','Location','southwest')
set(gcf, 'Color', 'w');

% Error figure
subplot(3,1,3)
semilogx(freq_Hz,error_a,'-o')
hold on
semilogx(freq_Hz,error_b,'-+')

semilogx(freq_Hz,error_c,'-^')

semilogx(freq_Hz,error_d,'-p')
semilogx(freq_Hz,error_e,'-*')
semilogx(freq_Hz,error_f,'-x')
semilogx(freq_Hz,error_g,'-s')
semilogx(freq_Hz,error_h,'-h')
semilogx(freq_Hz,error_j,'-d')
xlabel('Frequency (Hz)')
ylabel('Average Error Signal Amplitude')
title('Average Error Amplitudes')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','3nA input','Location','northwest')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','Location','northwest')
% legend('0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','Location','northwest')
legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','3nA input','4nA input','Location','northwest')

% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','Location','southwest')
set(gcf, 'Color', 'w');
%%
figure
subplot(3,1,1)
semilogx(freq_Hz,Kpsyngain_a,'-o')
hold on
semilogx(freq_Hz,Kpsyngain_b,'-+')
semilogx(freq_Hz,Kpsyngain_c,'-^')
semilogx(freq_Hz,Kpsyngain_d,'-p')
semilogx(freq_Hz,Kpsyngain_e,'-*')
semilogx(freq_Hz,Kpsyngain_f,'-x')
semilogx(freq_Hz,Kpsyngain_g,'-s')
semilogx(freq_Hz,Kpsyngain_h,'-h')
semilogx(freq_Hz,Kpsyngain_j,'-d')
% semilogx(freq_Hz,Kpsyngain_a+Kpstd_a,'--')
% semilogx(freq_Hz,Kpsyngain_j-Kpstd_j,'--')
% xlabel('Frequency (Hz)')
ylabel('Gain')
title('Kp Gain')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','3nA input','Location','southwest')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','Location','southwest')
set(gcf, 'Color', 'w');
ylim([1 2.5])

kdbase = 0.046;
subplot(3,1,2)
semilogx(freq_Hz,Kdsyngain_a*kdbase,'-o')
hold on
semilogx(freq_Hz,Kdsyngain_b*kdbase,'-+')
semilogx(freq_Hz,Kdsyngain_c*kdbase,'-^')
semilogx(freq_Hz,Kdsyngain_d*kdbase,'-p')
semilogx(freq_Hz,Kdsyngain_e*kdbase,'-*')
semilogx(freq_Hz,Kdsyngain_f*kdbase,'-x')
semilogx(freq_Hz,Kdsyngain_g*kdbase,'-s')
semilogx(freq_Hz,Kdsyngain_h*kdbase,'-h')
semilogx(freq_Hz,Kdsyngain_j*kdbase,'-d')
% semilogx(freq_Hz,Kdsyngain_a*kdbase+Kdstd_a*kdbase,'--')
% semilogx(freq_Hz,Kdsyngain_a*kdbase-Kdstd_a*kdbase,'--')
% xlabel('Frequency (Hz)')
ylabel('Gain')
title('Kd Gain')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','3nA input','Location','southwest')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','Location','southwest')
set(gcf, 'Color', 'w');

% Error figure
subplot(3,1,3)
semilogx(freq_Hz(1:5),error_a(1:5),'-o')
hold on
semilogx(freq_Hz,error_b,'-+')
semilogx(freq_Hz,error_c,'-^')
semilogx(freq_Hz,error_d,'-p')
semilogx(freq_Hz,error_e,'-*')
semilogx(freq_Hz,error_f,'-x')
semilogx(freq_Hz,error_g,'-s')
semilogx(freq_Hz,error_h,'-h')
semilogx(freq_Hz,error_j,'-d')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Average Error Signal Amplitude')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','3nA input','Location','northwest')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','Location','northwest')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','Location','northwest')
legend('1.5nA input','2nA input','3nA input','4nA input','Location','northwest')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','3nA input','4nA input','Location','northwest')
set(gcf, 'Color', 'w');
ylim([0 12])

%% extra plot for report
figure
subplot(3,1,1)
% semilogx(f_fftWNCL,magSys,'o')
% semilogx(f_fftWNCL,magWN,'Linewidth',1.4)
semilogx(freq_Hz,magRatioCL_a,'-o','Linewidth',1.2)
hold on
semilogx(freq_Hz,magRatioCL_b,'-+','Linewidth',1.2)
semilogx(freq_Hz,magRatioCL_c,'-^','Linewidth',1.2)
semilogx(freq_Hz,magRatioCL_d,'-p','Linewidth',1.2)
semilogx(freq_Hz,magRatioCL_e,'-*','Linewidth',1.2)
semilogx(freq_Hz,magRatioCL_f,'-x','Linewidth',1.2)
semilogx(freq_Hz,magRatioCL_g,'-s','Linewidth',1.2)
semilogx(freq_Hz,magRatioCL_h,'-h','Linewidth',1.2)
semilogx(freq_Hz,magRatioCL_j,'-d','Linewidth',1.2)
% semilogx(f_Human,magHuman,'Linewidth',1.2)
% semilogx(freq_Hz,magHigh,'--','Linewidth',1.4)
% semilogx(freq_Hz,magLow,'--','Linewidth',1.4)
xlim([0.1 2.25])
% fill(freq_Hz_flip,magflip,'b','FaceAlpha',0.1)
% xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('(a) Closed Loop System with Varied Auto Bode Inputs')

kdbase = 0.046;
subplot(3,1,2)
semilogx(freq_Hz,Kdsyngain_a*kdbase,'-o')
hold on
semilogx(freq_Hz,Kdsyngain_b*kdbase,'-+')
semilogx(freq_Hz,Kdsyngain_c*kdbase,'-^')
semilogx(freq_Hz,Kdsyngain_d*kdbase,'-p')
semilogx(freq_Hz,Kdsyngain_e*kdbase,'-*')
semilogx(freq_Hz,Kdsyngain_f*kdbase,'-x')
semilogx(freq_Hz,Kdsyngain_g*kdbase,'-s')
semilogx(freq_Hz,Kdsyngain_h*kdbase,'-h')
semilogx(freq_Hz,Kdsyngain_j*kdbase,'-d')
% semilogx(freq_Hz,Kdsyngain_a*kdbase+Kdstd_a*kdbase,'--')
% semilogx(freq_Hz,Kdsyngain_a*kdbase-Kdstd_a*kdbase,'--')
% xlabel('Frequency (Hz)')
ylabel('Gain')
title('(b) Kd Gain')
xlim([0.1 2.25])
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','3nA input','Location','southwest')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','Location','southwest')
set(gcf, 'Color', 'w');

% Error figure
subplot(3,1,3)
semilogx(freq_Hz(1:5),error_a(1:5),'-o')
hold on
semilogx(freq_Hz,error_b,'-+')
semilogx(freq_Hz,error_c,'-^')
semilogx(freq_Hz,error_d,'-p')
semilogx(freq_Hz,error_e,'-*')
semilogx(freq_Hz,error_f,'-x')
semilogx(freq_Hz,error_g,'-s')
semilogx(freq_Hz,error_h,'-h')
semilogx(freq_Hz,error_j,'-d')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('(c) Average Error Signal Amplitude')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','3nA input','Location','northwest')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','Location','northwest')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','Location','northwest')
% legend('1.5nA input','2nA input','3nA input','4nA input','Location','northwest')
legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','2nA input','3nA input','4nA input','Location','northwest')
set(gcf, 'Color', 'w');
ylim([0 12])
xlim([0.1 2.25])