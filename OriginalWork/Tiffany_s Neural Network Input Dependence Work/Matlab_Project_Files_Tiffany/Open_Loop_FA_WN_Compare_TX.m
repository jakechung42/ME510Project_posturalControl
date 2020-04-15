% Compare bodes for closed loop Controller overall, Kp and Kd circuits for varying
% inputs

clear
clc

% % Import Data
% % assumes all have the same frequency scale
A = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\PD_ABOL_12(0.1).dat');   % sys circuit outputs
freq = A(:,1);      % rad/s
magRatioOL_a = A(:,2);  
phaseOL_a = A(:,3);
OLstd_a = A(:,4);

B = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kp_ABOL_12(0.1).dat');   % Kp circuit outputs
magRatioKp_a = B(:,2);  
phaseKp_a = B(:,3);
error_a = B(:,4);
Kpstd_a = B(:,8);

C = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kd_ABOL_12(0.1).dat');   % Kd circuit outputs
magRatioKd_a = C(:,2);
phaseKd_a = C(:,3);

F1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kdsyngain_ABOL_12(0.1).dat');   % Kd circuit outputs
Kdsyngain_a = F1(:,2);

B2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\kdbase_ABOL_12(0.1).dat');   % Kd circuit outputs
kdbase_a = B2(:,2);

P1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kpsyngain_ABOL_12(0.1).dat');   % Kp syn gain
Kpsyngain_a = P1(:,2);
Kdsynstd_a = P1(:,4);

D = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\PD_ABOL_12(0.8).dat');   % sys circuit outputs
magRatioOL_b = D(:,2);  
phaseOL_b = D(:,3);
OLstd_b = D(:,4);

E = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kp_ABOL_12(0.8).dat');   % Kp circuit outputs
magRatioKp_b = E(:,2);  
phaseKp_b = E(:,3);
error_b = E(:,4);

F = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kd_ABOL_12(0.8).dat');   % Kd circuit outputs
magRatioKd_b = F(:,2);
phaseKd_b = F(:,3);
Kpstd_b = F(:,8);

G1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kdsyngain_ABOL_12(0.8).dat');   % Kd circuit outputs
Kdsyngain_b = G1(:,2);
Kdsynstd_b = G1(:,4);

C2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\kdbase_ABOL_12(0.8).dat');   % Kd circuit outputs
kdbase_b = C2(:,2);

Q1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kpsyngain_ABOL_12(0.8).dat');   % Kp syn gain
Kpsyngain_b = Q1(:,2);

G = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\PD_ABOL_12(4).dat');   % sys circuit outputs
magRatioOL_c = G(:,2);  
phaseOL_c = G(:,3);
OLstd_c = G(:,4);

H = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kp_ABOL_12(4).dat');   % Kp circuit outputs
magRatioKp_c = H(:,2);  
phaseKp_c = H(:,3);
error_c = H(:,4);

J = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kd_ABOL_12(4).dat');   % Kd circuit outputs
magRatioKd_c = J(:,2);
phaseKd_c = J(:,3);
Kpstd_c = J(:,8);

H1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kdsyngain_ABOL_12(4).dat');   % Kd circuit outputs
Kdsyngain_c = H1(:,2);
Kdsynstd_c = H1(:,4);

D2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\kdbase_ABOL_12(4).dat');   % Kd circuit outputs
kdbase_c = D2(:,2);

R1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kpsyngain_ABOL_12(4).dat');   % Kp syn gain
Kpsyngain_c = R1(:,2);

% K = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\PD_ABOL_4(0.75).dat');   % sys circuit outputs
% magRatioOL_d = K(:,2);  
% phaseOL_d = K(:,3);
% 
% L = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kp_ABOL_4(0.75).dat');   % Kp circuit outputs
% magRatioKp_d = L(:,2);  
% phaseKp_d = L(:,3);
% error_d = L(:,4);
% 
% P = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kd_ABOL_4(0.75).dat');   % Kd circuit outputs
% magRatioKd_d = P(:,2);
% phaseKd_d = P(:,3);
% 
% J1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kdsyngain_ABOL_4(0.75).dat');   % Kd circuit outputs
% Kdsyngain_d = J1(:,2);
% 
% E2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\kdbase_ABOL_4(0.75).dat');   % Kd circuit outputs
% kdbase_d = E2(:,2);
% 
% S1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kpsyngain_ABOL_4(0.75).dat');   % Kp syn gain
% Kpsyngain_d = S1(:,2);
% 
% Q = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\PD_ABOL_4(1).dat');   % sys circuit outputs
% magRatioOL_e = Q(:,2);  
% phaseOL_e = Q(:,3);
% 
% R = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kp_ABOL_4(1).dat');   % Kp circuit outputs
% magRatioKp_e = R(:,2);  
% phaseKp_e = R(:,3);
% error_e = R(:,4);
% 
% T = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kd_ABOL_4(1).dat');   % Kd circuit outputs
% magRatioKd_e = T(:,2);
% phaseKd_e = T(:,3);
% 
% K1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kdsyngain_ABOL_4(1).dat');   % Kd circuit outputs
% Kdsyngain_e = K1(:,2);
% 
% F2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\kdbase_ABOL_4(1).dat');   % Kd circuit outputs
% kdbase_e = F2(:,2);
% 
% T1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kpsyngain_ABOL_4(1).dat');   % Kp syn gain
% Kpsyngain_e = T1(:,2);
% 
% U = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\PD_ABOL_4(1.5).dat');   % sys circuit outputs
% magRatioOL_f = U(:,2);  
% phaseOL_f = U(:,3);
% 
% V = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kp_ABOL_4(1.5).dat');   % Kp circuit outputs
% magRatioKp_f = V(:,2);  
% phaseKp_f = V(:,3);
% error_f = V(:,4);
% 
% W = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kd_ABOL_4(1.5).dat');   % Kd circuit outputs
% magRatioKd_f = W(:,2);
% phaseKd_f = W(:,3);
% 
% L1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kdsyngain_ABOL_4(1.5).dat');   % Kd circuit outputs
% Kdsyngain_f = L1(:,2);
% 
% G2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\kdbase_ABOL_4(1.5).dat');   % Kd circuit outputs
% kdbase_f = G2(:,2);
% 
% U1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kpsyngain_ABOL_4(1.5).dat');   % Kp syn gain
% Kpsyngain_f = U1(:,2);
% 
% X = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\PD_ABOL_4(2).dat');   % sys circuit outputs
% magRatioOL_g = X(:,2);  
% phaseOL_g = X(:,3);
% 
% Y = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kp_ABOL_4(2).dat');   % Kp circuit outputs
% magRatioKp_g = Y(:,2);  
% phaseKp_g = Y(:,3);
% error_g = Y(:,4);
% 
% Z = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kd_ABOL_4(2).dat');   % Kd circuit outputs
% magRatioKd_g = Z(:,2);
% phaseKd_g = Z(:,3);
% 
% M1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kdsyngain_ABOL_4(2).dat');   % Kd circuit outputs
% Kdsyngain_g = M1(:,2);
% 
% H2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\kdbase_ABOL_4(2).dat');   % Kd circuit outputs
% kdbase_g = H2(:,2);
% 
% V1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kpsyngain_ABOL_4(2).dat');   % Kp syn gain
% Kpsyngain_g = V1(:,2);
% 
% A1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\PD_ABOL_4(3).dat');   % sys circuit outputs
% magRatioOL_h = A1(:,2);  
% phaseOL_h = A1(:,3);
% 
% B1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kp_ABOL_4(3).dat');   % Kp circuit outputs
% magRatioKp_h = B1(:,2);  
% phaseKp_h = B1(:,3);
% error_h = B1(:,4);
% 
% C1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kd_ABOL_4(3).dat');   % Kd circuit outputs
% magRatioKd_h = C1(:,2);
% phaseKd_h = C1(:,3);
% 
% N1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kdsyngain_ABOL_4(3).dat');   % Kd circuit outputs
% Kdsyngain_h = N1(:,2);
% 
% J2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\kdbase_ABOL_4(3).dat');   % Kd circuit outputs
% kdbase_h = J2(:,2);
% 
% W1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kpsyngain_ABOL_4(3).dat');   % Kp syn gain
% Kpsyngain_h = W1(:,2);
% 
% K2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\PD_ABOL_4(4).dat');   % sys circuit outputs
% magRatioOL_j = K2(:,2);  
% phaseOL_j = K2(:,3);
% 
% L2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kp_ABOL_4(4).dat');   % Kp circuit outputs
% magRatioKp_j = L2(:,2);  
% phaseKp_j = L2(:,3);
% error_j = L2(:,4);
% 
% M2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kd_ABOL_4(4).dat');   % Kd circuit outputs
% magRatioKd_j = M2(:,2);
% phaseKd_j = M2(:,3);
% 
% N2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kdsyngain_ABOL_4(4).dat');   % Kd circuit outputs
% Kdsyngain_j = N2(:,2);
% 
% P2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\kdbase_ABOL_4(4).dat');   % Kd circuit outputs
% kdbase_j = P2(:,2);
% 
% Q2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\ABOL_Outputs\Kpsyngain_ABOL_4(4).dat');   % Kp syn gain
% Kpsyngain_j = Q2(:,2);

% Closed Loop White Noise Data
% PD files
% Run a
D1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\ErrorN3_WNOL_64(0.2).dat');   % Desired input to NN
timeWNOL_a = D1(:,1);      % s
Error_N3WNOL_a = D1(:,2);  

E1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\NNOutput_WNOL_64(0.2).dat');   % Actual output (input to NN)
NNOutputWNOL_a = E1(:,2);  

% Kp files
X1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\ErrorScope_WNOL_64(0.2).dat');   % Desired input to NN
ErrorScopeWNOL_a = X1(:,2);  

Y1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\KpCircuitScope_WNOL_64(0.2).dat');   % Actual output (input to NN)
KpCircuitScopeWNOL_a = Y1(:,2); 

% Kd files
Z1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\kdbasegain_WNOL_64(0.2).dat');   % Desired input to NN
kdbasegainWNOL_a = Z1(:,2);  

A2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\KdCircuitScope_WNOL_64(0.2).dat');   % Actual output (input to NN)
KdCircuitScopeWNOL_a = A2(:,2); 

% Run b
A3 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\ErrorN3_WNOL_66(1.5).dat');   % Desired input to NN
timeWNOL_b = A3(:,1);      % s
Error_N3WNOL_b = A3(:,2);  

B3 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\NNOutput_WNOL_66(1.5).dat');   % Actual output (input to NN)
NNOutputWNOL_b = B3(:,2);  

% Kp files
C3 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\ErrorScope_WNOL_66(1.5).dat');   % Desired input to NN
ErrorScopeWNOL_b = C3(:,2);  

D3 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\KpCircuitScope_WNOL_66(1.5).dat');   % Actual output (input to NN)
KpCircuitScopeWNOL_b = D3(:,2); 

% Kd files
E3 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\kdbasegain_WNOL_66(1.5).dat');   % Desired input to NN
kdbasegainWNOL_b = E3(:,2);  

F3 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNOL_NNOutputs\KdCircuitScope_WNOL_66(1.5).dat');   % Actual output (input to NN)
KdCircuitScopeWNOL_b = F3(:,2); 

freq_Hz = freq/(2*pi);
% freq_Hz_flip = [freq_Hz,fliplr(freq_Hz)];

%% Open Loop Auto Bode Fits
s = tf('s');

fb1 = 5;         % cutoff frequency in Hz (base 2)
lowFiltDerivative = (1/(s/(2*pi*fb1)+1));  % low pass filter

fb2 = 80;         % cutoff frequency in Hz (base 80)
lowFiltNeuron = (1/(s/(2*pi*fb2)+1));  % low pass filter

% PD Circuit
tau = 0.004;
[num,den] = pade(tau,1);    % pade delay approximation first order (in order to use impulse command)
delay = tf(num,den);

Kpstart = 0.5; Kpend = 2.5; Kpdivs = 60;
Kdstart = 0.1; Kdend = 0.7; Kddivs = 60;
Kdbasestart = 0.01; Kdbaseend = 0.1;
% [KpPDFitLow,KdPDFitLow] = LSRfitOLPD_MagPhase_2Param(Kpstart,Kpend,Kpdivs,Kdstart,Kdend,Kddivs,freq_Hz,tau,magRatioOL_a,phaseOL_a)
% [KpPDFitMid,KdPDFitMid] = LSRfitOLPD_MagPhase_2Param(Kpstart,Kpend,Kpdivs,Kdstart,Kdend,Kddivs,freq_Hz,tau,magRatioOL_e,phaseOL_e)
% [KpPDFitHigh,KdPDFitHigh] = LSRfitOLPD_MagPhase_2Param(Kpstart,Kpend,Kpdivs,Kdstart,Kdend,Kddivs,freq_Hz,tau,magRatioOL_j,phaseOL_j)
% 
% GcHigh = lowFiltNeuron^4*(KpPDFitHigh + KdPDFitHigh*s*lowFiltDerivative);
% oltfHigh = delay*GcHigh;
% 
% [magHigh,phaseHigh] = bode(oltfHigh,freq);    % requires frequency in rad/s
% magHigh = squeeze(magHigh);
% phaseHigh = squeeze(phaseHigh);
% 
% GcMid = lowFiltNeuron^4*(KpPDFitMid + KdPDFitMid*s*lowFiltDerivative);
% oltfMid = delay*GcMid;
% 
% [magMid,phaseMid] = bode(oltfMid,freq);    % requires frequency in rad/s
% magMid = squeeze(magMid);
% phaseMid = squeeze(phaseMid);
% 
% GcLow = lowFiltNeuron^4*(KpPDFitLow + KdPDFitLow*s*lowFiltDerivative);
% oltfLow = delay*GcLow;
% 
% [magLow,phaseLow] = bode(oltfLow,freq);    % requires frequency in rad/s
% magLow = squeeze(magLow);
% phaseLow = squeeze(phaseLow);
% 
% % Kd Circuit
% [KdCircFitLow] = LSRfitKd_MagPhase_1Param(Kdstart,Kdend,Kddivs,freq_Hz,magRatioKd_a,phaseKd_a)
% [KdCircFitHigh] = LSRfitKd_MagPhase_1Param(Kdstart,Kdend,Kddivs,freq_Hz,magRatioKd_h,phaseKd_h)
% 
% oltfKdCircLow = lowFiltNeuron^2*KdCircFitLow*s*lowFiltDerivative;
% 
% [magKdLow,phaseKdLow] = bode(oltfKdCircLow,freq);    % requires frequency in rad/s
% magKdLow = squeeze(magKdLow);
% phaseKdLow = squeeze(phaseKdLow);
% 
% oltfKdCircHigh = lowFiltNeuron^2*KdCircFitHigh*s*lowFiltDerivative;
% 
% [magKdHigh,phaseKdHigh] = bode(oltfKdCircHigh,freq);    % requires frequency in rad/s
% magKdHigh = squeeze(magKdHigh);
% phaseKdHigh = squeeze(phaseKdHigh);
% 
% % kd base Fits
% [KdbaseFitLow] = LSRfitKd_MagPhase_1Param(Kdstart,Kdend,Kddivs,freq_Hz,magRatioKd_a,phaseKd_a)
% [KdbaseFitHigh] = LSRfitKd_MagPhase_1Param(Kdstart,Kdend,Kddivs,freq_Hz,magRatioKd_h,phaseKd_h)
% 
% oltfKdbaseLow = lowFiltNeuron^2*KdbaseFitLow*s*lowFiltDerivative;
% 
% [magKdbaseLow,phaseKdbaseLow] = bode(oltfKdbaseLow,freq);    % requires frequency in rad/s
% magKdbaseLow = squeeze(magKdbaseLow);
% 
% oltfKdbaseHigh = lowFiltNeuron^2*KdbaseFitHigh*s*lowFiltDerivative;
% 
% [magKdbaseHigh,phaseKdbaseHigh] = bode(oltfKdbaseHigh,freq);    % requires frequency in rad/s
% magKdbaseHigh = squeeze(magKdbaseHigh);
% 
% % Kp Circuit
% [KpCircFitLow] = LSRfitKp_MagPhase_1Param(Kpstart,Kpend,Kpdivs,freq_Hz,magRatioKp_a,phaseKp_a)
% [KpCircFitHigh] = LSRfitKp_MagPhase_1Param(Kpstart,Kpend,Kpdivs,freq_Hz,magRatioKp_h,phaseKp_h)
% 
% oltfKpCircLow = lowFiltNeuron^2*KpCircFitLow;
% 
% [magKpLow,phaseKpLow] = bode(oltfKpCircLow,freq);    % requires frequency in rad/s
% magKpLow = squeeze(magKpLow);
% 
% oltfKpCircHigh = lowFiltNeuron^2*KpCircFitHigh;
% 
% [magKpHigh,phaseKpHigh] = bode(oltfKpCircHigh,freq);    % requires frequency in rad/s
% magKpHigh = squeeze(magKpHigh);

%% Process White Noise Data
% Run a
dtWNOL_a = timeWNOL_a(3,1)-timeWNOL_a(2,1);
TsWNOL_a = 1/dtWNOL_a;
numSegmentsWNOL = 10;   % # segments to divide signal into for fft averaging
movAvgWinTimeWNOL = 0.5;   % window size for moving average in seconds

[f_fftWNOL_a,magSys_a,phaseSys_a,AmpinSys_a,AmpoutSys_a] = bodebyFFT(Error_N3WNOL_a,NNOutputWNOL_a,TsWNOL_a,numSegmentsWNOL,movAvgWinTimeWNOL);
[f_fftWNOL_a,magKpCirc_a,phaseKpCirc_a,AmpinKpCirc_a,AmpoutKpCirc_a] = bodebyFFT(ErrorScopeWNOL_a,KpCircuitScopeWNOL_a,TsWNOL_a,numSegmentsWNOL,movAvgWinTimeWNOL);
[f_fftWNOL_a,magKdCirc_a,phaseKdCirc_a,AmpinKdCirc_a,AmpoutKdCirc_a] = bodebyFFT(ErrorScopeWNOL_a,KdCircuitScopeWNOL_a,TsWNOL_a,numSegmentsWNOL,movAvgWinTimeWNOL);
[f_fftWNOL_a,magkdbase_a,phasekdbase_a,Ampinkdbasegain_a,Ampoutkdbasegain_a] = bodebyFFT(ErrorScopeWNOL_a,kdbasegainWNOL_a,TsWNOL_a,numSegmentsWNOL,movAvgWinTimeWNOL);

[KpPDFitWN_a,KdPDFitWN_a,PDstdev_a] = LSRfitOLPD_MagPhase_2Param_Error(Kpstart,Kpend,Kpdivs,Kdstart,Kdend,Kddivs,f_fftWNOL_a(1:floor(end/4)),tau,magSys_a(1:floor(end/4)),phaseSys_a(1:floor(end/4)));
[KpCircFitWN_a,Kpstdev_a] = LSRfitKp_MagPhase_1Param_Error(Kpstart,Kpend,Kpdivs,f_fftWNOL_a(1:floor(end/4)),magKpCirc_a(1:floor(end/4)),phaseKpCirc_a(1:floor(end/4)));
[KdCircFitWN_a,Kdstdev_a] = LSRfitKd_MagPhase_1Param_Error(Kdstart,Kdend,Kddivs,f_fftWNOL_a(1:floor(end/4)),magKdCirc_a(1:floor(end/4)),phaseKdCirc_a(1:floor(end/4)));
[KdbaseFitWN_a,Kdbasestdev_a] = LSRfitKd_MagPhase_1Param_Error(Kdbasestart,Kdbaseend,Kddivs,f_fftWNOL_a(1:floor(end/4)),magkdbase_a(1:floor(end/4)),phasekdbase_a(1:floor(end/4)));

% Run b
dtWNOL_b = timeWNOL_b(3,1)-timeWNOL_b(2,1);
TsWNOL_b = 1/dtWNOL_b;
numSegmentsWNOL = 10;   % # segments to divide signal into for fft averaging
movAvgWinTimeWNOL = 0.5;   % window size for moving average in seconds

[f_fftWNOL_b,magSys_b,phaseSys_b,AmpinSys_b,AmpoutSys_b] = bodebyFFT(Error_N3WNOL_b,NNOutputWNOL_b,TsWNOL_b,numSegmentsWNOL,movAvgWinTimeWNOL);
[f_fftWNOL_b,magKpCirc_b,phaseKpCirc_b,AmpinKpCirc_b,AmpoutKpCirc_b] = bodebyFFT(ErrorScopeWNOL_b,KpCircuitScopeWNOL_b,TsWNOL_b,numSegmentsWNOL,movAvgWinTimeWNOL);
[f_fftWNOL_b,magKdCirc_b,phaseKdCirc_b,AmpinKdCirc_b,AmpoutKdCirc_b] = bodebyFFT(ErrorScopeWNOL_b,KdCircuitScopeWNOL_b,TsWNOL_b,numSegmentsWNOL,movAvgWinTimeWNOL);
[f_fftWNOL_b,magkdbase_b,phasekdbase_b,Ampinkdbasegain_b,Ampoutkdbasegain_b] = bodebyFFT(ErrorScopeWNOL_b,kdbasegainWNOL_b,TsWNOL_b,numSegmentsWNOL,movAvgWinTimeWNOL);

[KpPDFitWN_b,KdPDFitWN_b,PDstdev_b] = LSRfitOLPD_MagPhase_2Param_Error(Kpstart,Kpend,Kpdivs,Kdstart,Kdend,Kddivs,f_fftWNOL_b(1:floor(end/4)),tau,magSys_b(1:floor(end/4)),phaseSys_b(1:floor(end/4)));
[KpCircFitWN_b,Kpstdev_b] = LSRfitKp_MagPhase_1Param_Error(Kpstart,Kpend,Kpdivs,f_fftWNOL_b(1:floor(end/4)),magKpCirc_b(1:floor(end/4)),phaseKpCirc_b(1:floor(end/4)));
[KdCircFitWN_b,Kdstdev_b] = LSRfitKd_MagPhase_1Param_Error(Kdstart,Kdend,Kddivs,f_fftWNOL_b(1:floor(end/4)),magKdCirc_b(1:floor(end/4)),phaseKdCirc_b(1:floor(end/4)));
[KdbaseFitWN_b,Kdbasestdev_b] = LSRfitKd_MagPhase_1Param_Error(Kdbasestart,Kdbaseend,Kddivs,f_fftWNOL_b(1:floor(end/4)),magkdbase_b(1:floor(end/4)),phasekdbase_b(1:floor(end/4)));

%%
% magflip = [magHigh,fliplr(magLow)];
% 
% patch(freq_Hz_flip,magflip,'b','FaceAlpha',0.1)
% figure
% hold on
% area(freq_Hz,magLow,'FaceColor','b')
% area(freq_Hz,magHigh,'FaceColor','w')
%%
% Open loop fit to whitenoise data
% Run a
GcKpWN_a = lowFiltNeuron^2*KpCircFitWN_a;

GcKdWN_a = lowFiltNeuron^2*KdCircFitWN_a*s*lowFiltDerivative;

GckdbaseWN_a = lowFiltNeuron^2*KdbaseFitWN_a*s*lowFiltDerivative;

GcWN_a = lowFiltNeuron^4*(KpPDFitWN_a + KdPDFitWN_a*s*lowFiltDerivative);
oltfWN_a = delay*GcWN_a;

GcKpKdWN_a = lowFiltNeuron^4*(KpCircFitWN_a + KdCircFitWN_a*s*lowFiltDerivative);
oltfKpKdWN_a = delay*GcKpKdWN_a;

f_fftWNCL_rad_a = f_fftWNOL_a*2*pi;
[magKpWN_a,phaseKpWN_a] = bode(GcKpWN_a,f_fftWNCL_rad_a);    % requires frequency in rad/s
magKpWN_a = squeeze(magKpWN_a);

[magKdWN_a,phaseKdWN_a] = bode(GcKdWN_a,f_fftWNCL_rad_a);    % requires frequency in rad/s
magKdWN_a = squeeze(magKdWN_a);

[magKdbaseWN_a,phaseKdbaseWN_a] = bode(GckdbaseWN_a,f_fftWNCL_rad_a);    % requires frequency in rad/s
magKdbaseWN_a = squeeze(magKdbaseWN_a);

[magWN_a,phaseWN_a] = bode(oltfWN_a,f_fftWNCL_rad_a);    % requires frequency in rad/s
magWN_a = squeeze(magWN_a);
phaseWN_a = squeeze(phaseWN_a);

[magKpKdWN_a,phaseKpKdWN_a] = bode(oltfKpKdWN_a,f_fftWNCL_rad_a);    % requires frequency in rad/s
magKpKdWN_a = squeeze(magKpKdWN_a);

% Run b
GcKpWN_b = lowFiltNeuron^2*KpCircFitWN_b;

GcKdWN_b = lowFiltNeuron^2*KdCircFitWN_b*s*lowFiltDerivative;

GckdbaseWN_b = lowFiltNeuron^2*KdbaseFitWN_b*s*lowFiltDerivative;

GcWN_b = lowFiltNeuron^4*(KpPDFitWN_b + KdPDFitWN_b*s*lowFiltDerivative);
oltfWN_b = delay*GcWN_b;

GcKpKdWN_b = lowFiltNeuron^4*(KpCircFitWN_b + KdCircFitWN_b*s*lowFiltDerivative);
oltfKpKdWN_b = delay*GcKpKdWN_b;

f_fftWNCL_rad_b = f_fftWNOL_b*2*pi;
[magKpWN_b,phaseKpWN_b] = bode(GcKpWN_b,f_fftWNCL_rad_b);    % requires frequency in rad/s
magKpWN_b = squeeze(magKpWN_b);

[magKdWN_b,phaseKdWN_b] = bode(GcKdWN_b,f_fftWNCL_rad_b);    % requires frequency in rad/s
magKdWN_b = squeeze(magKdWN_b);

[magKdbaseWN_b,phaseKdbaseWN_b] = bode(GckdbaseWN_b,f_fftWNCL_rad_b);    % requires frequency in rad/s
magKdbaseWN_b = squeeze(magKdbaseWN_b);

[magWN_b,phaseWN_b] = bode(oltfWN_b,f_fftWNCL_rad_b);    % requires frequency in rad/s
magWN_b = squeeze(magWN_b);
phaseWN_b = squeeze(phaseWN_b);

[magKpKdWN_b,phaseKpKdWN_b] = bode(oltfKpKdWN_b,f_fftWNCL_rad_b);    % requires frequency in rad/s
magKpKdWN_b = squeeze(magKpKdWN_b);

%% Plots

% Run a
% Open Loop PD figure
figure
semilogx(f_fftWNOL_a,magSys_a,'o')
hold on
semilogx(f_fftWNOL_a,magWN_a,'Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_a,'-s','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_a+3*OLstd_a,'--','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_a-3*OLstd_a,'--','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_b,'-+','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_b+3*OLstd_b,'--','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_b-3*OLstd_b,'--','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_c,'-^','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_c+3*OLstd_c,'--','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_c-3*OLstd_c,'--','Linewidth',1.4)
% fill(freq_Hz_flip,magflip,'b','FaceAlpha',0.1)
xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('Open Loop System Transfer Function with Variable Inputs')
% legend('White Noise PD Data',['Fit Kp = ' num2str(KpPDFitWN_a) ',  Kd = ' num2str(KdPDFitWN_a)],'0.1nA input','0.1nA input+3std','0.1nA input-3std','0.8nA input','0.8nA input+3std','0.8nA input-3std','4nA input','4nA input+3std','4nA input-3std','Location','northwest')
legend('White Noise PD Data',['Fit Kp = ' num2str(KpPDFitWN_a) ',  Kd = ' num2str(KdPDFitWN_a) ',  PD std = ' num2str(PDstdev_a,2)],'Location','northwest')
% ylim([-40 20])
xlim([0.1 100])
set(gcf, 'Color', 'w');
% export_fig('C:\Users\Tiffany\Desktop\Sim_Outputs\Figs\OLABWN_MultiInput_KpbaseKd095_Error', '-bmp','-q101','-nocrop','-a1');

% kd base gain figure
figure
semilogx(f_fftWNOL_a,magkdbase_a,'o')
hold on
semilogx(f_fftWNOL_a,magKdbaseWN_a,'Linewidth',1.4)
% semilogx(freq_Hz,kdbase_a,'-o','Linewidth',1.4)
% semilogx(freq_Hz,kdbase_b,'-+','Linewidth',1.4)
% semilogx(freq_Hz,kdbase_c,'-^','Linewidth',1.4)
xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('Baseline kd Gain (by fitting TF) with Variable Inputs')
legend('Whitenoise kd Baseline Data',['Fit baseKd = ' num2str(KdbaseFitWN_a)],'Location','northwest')
xlim([0.1 100])
ylim([0 2.5])
set(gcf, 'Color', 'w');

% Kd figure
figure
subplot(2,1,1)
semilogx(f_fftWNOL_a,magKdCirc_a,'o')
hold on
semilogx(f_fftWNOL_a,magKdWN_a)
% semilogx(freq_Hz,magRatioKd_a,'-o')
% semilogx(freq_Hz,magRatioKd_b,'-+')
% semilogx(freq_Hz,magRatioKd_c,'-^')
xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('Kd Circuit Transfer Function Variable Inputs')
legend('White Noise Kd Data',['Fit Kd = ' num2str(KdCircFitWN_a)],'Location','northwest')
% ylim([0 5])
xlim([0.1 100])
set(gcf, 'Color', 'w');

% Kp figure
subplot(2,1,2)
semilogx(f_fftWNOL_a,magKpCirc_a,'o')
hold on
semilogx(f_fftWNOL_a,magKpWN_a)
% semilogx(freq_Hz,magRatioKp_a,'-o','Linewidth',1.4)
% semilogx(freq_Hz,magRatioKp_b,'-+','Linewidth',1.4)
% semilogx(freq_Hz,magRatioKp_c,'-^','Linewidth',1.4)
% semilogx(freq_Hz,magKp)
xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('Kp Circuit Transfer Function Variable Inputs')
legend('White Noise Kp Data',['Fit Kd = ' num2str(KpCircFitWN_a)],'Location','northwest')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','Location','southwest')
xlim([0.1 100])
set(gcf, 'Color', 'w');
%%
% Run b
% Open Loop PD figure
figure
semilogx(f_fftWNOL_b,magSys_b,'o')
hold on
% semilogx(freq_Hz,magRatioOL_a,'-s','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_b,'-+','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_c,'-^','Linewidth',1.4)
semilogx(f_fftWNOL_b,magWN_b,'Linewidth',1.4)
% fill(freq_Hz_flip,magflip,'b','FaceAlpha',0.1)
xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('Open Loop System Transfer Function with Variable Inputs')
legend('White Noise PD Data',['Fit Kp = ' num2str(KpPDFitWN_b) ',  Kd = ' num2str(KdPDFitWN_b)],'Location','northwest')
% ylim([-40 20])
xlim([0.1 100])
set(gcf, 'Color', 'w');
% export_fig('C:\Users\Tiffany\Desktop\Sim_Outputs\Figs\OLABWN_MultiInput_KpbaseKd095_Error', '-bmp','-q101','-nocrop','-a1');

% kd base gain figure
figure
semilogx(f_fftWNOL_b,magkdbase_b,'o')
hold on
% semilogx(freq_Hz,kdbase_a,'-o','Linewidth',1.4)
% semilogx(freq_Hz,kdbase_b,'-+','Linewidth',1.4)
% semilogx(freq_Hz,kdbase_c,'-^','Linewidth',1.4)
semilogx(f_fftWNOL_b,magKdbaseWN_b,'Linewidth',1.4)
xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('Baseline kd Gain (by fitting TF) with Variable Inputs')
legend('Whitenoise kd Baseline Data',['Fit baseKd = ' num2str(KdbaseFitWN_b)],'Location','northwest')
xlim([0.1 100])
ylim([0 2.5])
set(gcf, 'Color', 'w');

% Kd figure
figure
subplot(2,1,1)
semilogx(f_fftWNOL_b,magKdCirc_b,'o')
hold on
% semilogx(freq_Hz,magRatioKd_a,'-o')
% semilogx(freq_Hz,magRatioKd_b,'-+')
% semilogx(freq_Hz,magRatioKd_c,'-^')
semilogx(f_fftWNOL_b,magKdWN_b)
xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('Kd Circuit Transfer Function Variable Inputs')
legend('White Noise Kd Data',['Fit Kd = ' num2str(KdCircFitWN_b)],'Location','northwest')
% ylim([0 5])
xlim([0.1 100])
set(gcf, 'Color', 'w');

% Kp figure
subplot(2,1,2)
semilogx(f_fftWNOL_b,magKpCirc_b,'o')
hold on
% semilogx(freq_Hz,magRatioKp_a,'-o','Linewidth',1.4)
% semilogx(freq_Hz,magRatioKp_b,'-+','Linewidth',1.4)
% semilogx(freq_Hz,magRatioKp_c,'-^','Linewidth',1.4)
semilogx(f_fftWNOL_b,magKpWN_b)
% semilogx(freq_Hz,magKp)
xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('Kp Circuit Transfer Function Variable Inputs')
legend('White Noise Kp Data',['Fit Kd = ' num2str(KpCircFitWN_b)],'Location','northwest')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','Location','southwest')
xlim([0.1 100])
set(gcf, 'Color', 'w');
%%
% Combo Fig
% Open Loop PD figure
figure
% semilogx(f_fftWNOL_b,magSys_b,'o')
semilogx(f_fftWNOL_a,magWN_a,'Linewidth',1.4,'Color','k')
hold on
semilogx(f_fftWNOL_a,magWN_a+PDstdev_a,'-.','Color','k')
semilogx(f_fftWNOL_a,magWN_a-PDstdev_a,'-.','Color','k')
% semilogx(f_fftWNOL_a,magKpKdWN_a,'--','Linewidth',1.4,'Color','k')
semilogx(f_fftWNOL_b,magWN_b,'Linewidth',1.4,'Color','g')
semilogx(f_fftWNOL_b,magWN_b+PDstdev_b,'-.','Color','g')
semilogx(f_fftWNOL_b,magWN_b-PDstdev_b,'-.','Color','g')
% semilogx(f_fftWNOL_b,magKpKdWN_b,'--','Linewidth',1.4,'Color','g')
% semilogx(freq_Hz,magRatioOL_a,'-s','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_a+3*OLstd_a,'--','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_a-3*OLstd_a,'--','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_b,'-+','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_b+3*OLstd_b,'--','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_b-3*OLstd_b,'--','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_c,'-^','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_c+3*OLstd_c,'--','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_c-3*OLstd_c,'--','Linewidth',1.4)
% fill(freq_Hz_flip,magflip,'b','FaceAlpha',0.1)
xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('(a) Open Loop System: Varied Mean Inputs')
% legend(['Run a: Fit Kp = ' num2str(KpPDFitWN_a,2) ',  Kd = ' num2str(KdPDFitWN_a,2)],['Run a: +1std = ' num2str(PDstdev_a,2)],'Run a: -1std','Run a: Fit from Kp/Kd Circuits',['Run b: Fit Kp = ' num2str(KpPDFitWN_b,2) ',  Kd = ' num2str(KdPDFitWN_b,2)],['Run b: +1std = ' num2str(PDstdev_b,2)],'Run b: -1std','Run b: Fit from Kp/Kd Circuits','Location','northwest')
legend(['0.2 nA Mean Amp: Fit Kp = ' num2str(KpPDFitWN_a,2) ',  Kd = ' num2str(KdPDFitWN_a,2)],['0.2 nA: +1std = ' num2str(PDstdev_a,2)],'0.2 nA: -1std',['1.5 nA Mean Amp: Fit Kp = ' num2str(KpPDFitWN_b,2) ',  Kd = ' num2str(KdPDFitWN_b,2)],['1.5 nA: +1std = ' num2str(PDstdev_b,2)],'1.5 nA: -1std','Location','northwest')

% legend(['Run a: Fit Kp = ' num2str(KpPDFitWN_a,2) ',  Kd = ' num2str(KdPDFitWN_a,2)],['Run b: Fit Kp = ' num2str(KpPDFitWN_b,2) ',  Kd = ' num2str(KdPDFitWN_b,2)],'Run a: Fit from Kp/Kd Circuits','Run b: Fit from Kp/Kd Circuits','0.1nA input','0.1nA input+3std','0.1nA input-3std','0.8nA input','0.8nA input+3std','0.8nA input-3std','4nA input','4nA input+3std','4nA input-3std','Location','northwest')
ylim([0 12])
xlim([0.1 100])
set(gcf, 'Color', 'w');
% export_fig('C:\Users\Tiffany\Desktop\Sim_Outputs\Figs\OLABWN_MultiInput_KpbaseKd095_Error', '-bmp','-q101','-nocrop','-a1');

% kd base gain figure
figure
% semilogx(f_fftWNOL_b,magkdbase_b,'o')
semilogx(f_fftWNOL_a,magKdbaseWN_a,'Linewidth',1.4,'Color','k')
hold on
semilogx(f_fftWNOL_a,magKdbaseWN_a+Kdbasestdev_a,'-.','Color','k')
semilogx(f_fftWNOL_a,magKdbaseWN_a-Kdbasestdev_a,'-.','Color','k')
semilogx(f_fftWNOL_b,magKdbaseWN_b,'Linewidth',1.4,'Color','m')
semilogx(f_fftWNOL_b,magKdbaseWN_b+Kdbasestdev_b,'-.','Color','m')
semilogx(f_fftWNOL_b,magKdbaseWN_b-Kdbasestdev_b,'-.','Color','m')
% semilogx(freq_Hz,kdbase_a,'-o','Linewidth',1.4)
% semilogx(freq_Hz,kdbase_b,'-+','Linewidth',1.4)
% semilogx(freq_Hz,kdbase_c,'-^','Linewidth',1.4)
xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('Baseline kd Transfer Function')
legend(['0.2 nA: Fit baseKd = ' num2str(KdbaseFitWN_a,2)],['0.2 nA: +1std = ' num2str(Kdbasestdev_a,2)],'0.2 nA: -1std',['1.5 nA: Fit baseKd = ' num2str(KdbaseFitWN_b,2)],['1.5 nA: +1std = ' num2str(Kdbasestdev_b,2)],'1.5 nA: -1std','Location','northwest')
xlim([0.1 100])
% ylim([0 2.5])
set(gcf, 'Color', 'w');

% Kd figure
figure
subplot(2,1,1)
semilogx(f_fftWNOL_a,magKdWN_a,'Linewidth',1.4,'Color','k')
hold on
semilogx(f_fftWNOL_a,magKdWN_a+Kdstdev_a,'-.','Color','k')
semilogx(f_fftWNOL_a,magKdWN_a-Kdstdev_a,'-.','Color','k')
semilogx(f_fftWNOL_b,magKdWN_b,'Linewidth',1.4,'Color','g')
semilogx(f_fftWNOL_b,magKdWN_b+Kdstdev_b,'-.','Color','g')
semilogx(f_fftWNOL_b,magKdWN_b-Kdstdev_b,'-.','Color','g')
% semilogx(f_fftWNOL_b,magKdCirc_b,'o')
% semilogx(freq_Hz,magRatioKd_a,'-o')
% semilogx(freq_Hz,magRatioKd_b,'-+')
% semilogx(freq_Hz,magRatioKd_c,'-^')
% xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('(b) Kd Circuit')
legend(['0.2 nA: Fit Kd = ' num2str(KdCircFitWN_a,2)],['0.2 nA: +1std = ' num2str(PDstdev_a,2)],'0.2 nA: -1std',['1.5 nA: Fit Kd = ' num2str(KdCircFitWN_b,2)],['1.5 nA: +1std = ' num2str(PDstdev_b,2)],'1.5 nA: -1std','Location','northwest')
ylim([0 12])
xlim([0.1 100])
set(gcf, 'Color', 'w');

% Kp figure
subplot(2,1,2)
% semilogx(f_fftWNOL_b,magKpCirc_b,'o')
semilogx(f_fftWNOL_a,magKpWN_a,'Linewidth',1.4,'Color','k')
hold on
semilogx(f_fftWNOL_a,magKpWN_a+Kpstdev_a,'-.','Color','k')
semilogx(f_fftWNOL_a,magKpWN_a-Kpstdev_a,'-.','Color','k')
semilogx(f_fftWNOL_b,magKpWN_b,'Linewidth',1.4,'Color','g')
semilogx(f_fftWNOL_b,magKpWN_b+Kpstdev_b,'-.','Color','g')
semilogx(f_fftWNOL_b,magKpWN_b-Kpstdev_b,'-.','Color','g')
% semilogx(freq_Hz,magRatioKp_a,'-o','Linewidth',1.4)
% semilogx(freq_Hz,magRatioKp_b,'-+','Linewidth',1.4)
% semilogx(freq_Hz,magRatioKp_c,'-^','Linewidth',1.4)
% semilogx(freq_Hz,magKp)
xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('(c) Kp Circuit')
legend(['0.2 nA: Fit Kp = ' num2str(KpCircFitWN_a,2)],['0.2 nA: +1std = ' num2str(Kpstdev_a,2)],'0.2 nA: -1std',['1.5 nA: Fit Kp = ' num2str(KpCircFitWN_b,2)],['1.5 nA: +1std = ' num2str(Kpstdev_b,2)],'1.5 nA: -1std','Location','southwest')
% legend('0.1nA input','0.25nA input','0.5nA input','0.75nA input','1nA input','1.5nA input','Location','southwest')
xlim([0.1 100])
set(gcf, 'Color', 'w');

%% Auto Bode Figs
figure
% semilogx(f_fftWNOL_b,magSys_b,'o')
% semilogx(f_fftWNOL_a,magWN_a,'Linewidth',1.4)
% semilogx(f_fftWNOL_a,magKpKdWN_a,'--','Linewidth',1.4)
% semilogx(f_fftWNOL_b,magWN_b,'Linewidth',1.4)
% semilogx(f_fftWNOL_b,magWN_b+PDstdev_b,'-.')
% semilogx(f_fftWNOL_b,magWN_b-PDstdev_b,'-.')
% semilogx(f_fftWNOL_b,magKpKdWN_b,'--','Linewidth',1.4)
semilogx(freq_Hz,magRatioOL_a,'-s','Linewidth',1.4,'Color','b')
hold on
semilogx(freq_Hz,magRatioOL_a+3*OLstd_a,'--','Linewidth',1.4,'Color','b')
semilogx(freq_Hz,magRatioOL_a-3*OLstd_a,'--','Linewidth',1.4,'Color','b')
semilogx(freq_Hz,magRatioOL_b,'-+','Linewidth',1.4,'Color','g')
semilogx(freq_Hz,magRatioOL_b+3*OLstd_b,'--','Linewidth',1.4,'Color','g')
semilogx(freq_Hz,magRatioOL_b-3*OLstd_b,'--','Linewidth',1.4,'Color','g')
semilogx(freq_Hz,magRatioOL_c,'-^','Linewidth',1.4,'Color','k')
semilogx(freq_Hz,magRatioOL_c+3*OLstd_c,'--','Linewidth',1.4,'Color','k')
semilogx(freq_Hz,magRatioOL_c-3*OLstd_c,'--','Linewidth',1.4,'Color','k')
% fill(freq_Hz_flip,magflip,'b','FaceAlpha',0.1)
xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('(a) Open Loop System with Varied Auto Bode Inputs')
% legend(['Run a: Fit Kp = ' num2str(KpPDFitWN_a,2) ',  Kd = ' num2str(KdPDFitWN_a,2)],['Run a: +1std = ' num2str(PDstdev_a,2)],'Run a: -1std','Run a: Fit from Kp/Kd Circuits',['Run b: Fit Kp = ' num2str(KpPDFitWN_b,2) ',  Kd = ' num2str(KdPDFitWN_b,2)],['Run b: +1std = ' num2str(PDstdev_b,2)],'Run b: -1std','Run b: Fit from Kp/Kd Circuits','Location','northwest')
legend('0.1nA input','0.1nA input+3std','0.1nA input-3std','0.8nA input','0.8nA input+3std','0.8nA input-3std','4nA input','4nA input+3std','4nA input-3std','Location','northwest')
ylim([0 12])
xlim([0.5 100])
set(gcf, 'Color', 'w');
% export_fig('C:\Users\Tiffany\Desktop\Sim_Outputs\Figs\OLABWN_MultiInput_KpbaseKd095_Error', '-bmp','-q101','-nocrop','-a1');

% Kd figure
figure
% semilogx(f_fftWNOL_a,magKdWN_a,'Linewidth',1.4)
% semilogx(f_fftWNOL_a,magKdWN_a+Kdstdev_a,'-.')
% semilogx(f_fftWNOL_a,magKdWN_a-Kdstdev_a,'-.')
% semilogx(f_fftWNOL_b,magKdWN_b,'Linewidth',1.4)
% semilogx(f_fftWNOL_b,magKdWN_b+Kdstdev_b,'-.')
% semilogx(f_fftWNOL_b,magKdWN_b-Kdstdev_b,'-.')
% semilogx(f_fftWNOL_b,magKdCirc_b,'o')
% semilogx(freq_Hz,Kdsyngain_a,'-o')
semilogx(freq_Hz,Kdsyngain_a,'-s','Linewidth',1.4,'Color','b')
hold on
% semilogx(freq_Hz,Kdsyngain_b,'-+')
% semilogx(freq_Hz,Kdsyngain_c,'-^')
semilogx(freq_Hz,Kdsyngain_a+3*Kdsynstd_a,'--','Linewidth',1.4,'Color','b')
semilogx(freq_Hz,Kdsyngain_a-3*Kdsynstd_a,'--','Linewidth',1.4,'Color','b')
semilogx(freq_Hz,Kdsyngain_b,'-+','Linewidth',1.4,'Color','g')
semilogx(freq_Hz,Kdsyngain_b+3*Kdsynstd_b,'--','Linewidth',1.4,'Color','g')
semilogx(freq_Hz,Kdsyngain_b-3*Kdsynstd_b,'--','Linewidth',1.4,'Color','g')
semilogx(freq_Hz,Kdsyngain_c,'-^','Linewidth',1.4,'Color','k')
semilogx(freq_Hz,Kdsyngain_c+3*Kdsynstd_c,'--','Linewidth',1.4,'Color','k')
semilogx(freq_Hz,Kdsyngain_c-3*Kdsynstd_c,'--','Linewidth',1.4,'Color','k')
xlabel('Frequency (Hz)')
ylabel('Gain (abs)')
title('(b) Kd Synaptic Gain')
legend('0.1nA input','0.1nA input+3std','0.1nA input-3std','0.8nA input','0.8nA input+3std','0.8nA input-3std','4nA input','4nA input+3std','4nA input-3std','Location','southwest')

ylim([4 7])
xlim([0.1 100])
set(gcf, 'Color', 'w');

%% Auto Bode Figs (1 std)
figure
% semilogx(f_fftWNOL_b,magSys_b,'o')
% semilogx(f_fftWNOL_a,magWN_a,'Linewidth',1.4)
% semilogx(f_fftWNOL_a,magKpKdWN_a,'--','Linewidth',1.4)
% semilogx(f_fftWNOL_b,magWN_b,'Linewidth',1.4)
% semilogx(f_fftWNOL_b,magWN_b+PDstdev_b,'-.')
% semilogx(f_fftWNOL_b,magWN_b-PDstdev_b,'-.')
% semilogx(f_fftWNOL_b,magKpKdWN_b,'--','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_a,'-s','Linewidth',1.4,'Color','b')

% semilogx(freq_Hz,magRatioOL_a+3*OLstd_a,'--','Linewidth',1.4,'Color','b')
% semilogx(freq_Hz,magRatioOL_a-3*OLstd_a,'--','Linewidth',1.4,'Color','b')
semilogx(freq_Hz,magRatioOL_b,'-+','Linewidth',1.4,'Color','k')
hold on
semilogx(freq_Hz,magRatioOL_b+OLstd_b,'--','Linewidth',1.4,'Color','k')
semilogx(freq_Hz,magRatioOL_b-OLstd_b,'--','Linewidth',1.4,'Color','k')
% semilogx(freq_Hz,magRatioOL_c,'-^','Linewidth',1.4,'Color','k')
% semilogx(freq_Hz,magRatioOL_c+3*OLstd_c,'--','Linewidth',1.4,'Color','k')
% semilogx(freq_Hz,magRatioOL_c-3*OLstd_c,'--','Linewidth',1.4,'Color','k')
% fill(freq_Hz_flip,magflip,'b','FaceAlpha',0.1)
xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('Auto Bode Method Open Loop + 1 standard deviation')
% legend(['Run a: Fit Kp = ' num2str(KpPDFitWN_a,2) ',  Kd = ' num2str(KdPDFitWN_a,2)],['Run a: +1std = ' num2str(PDstdev_a,2)],'Run a: -1std','Run a: Fit from Kp/Kd Circuits',['Run b: Fit Kp = ' num2str(KpPDFitWN_b,2) ',  Kd = ' num2str(KdPDFitWN_b,2)],['Run b: +1std = ' num2str(PDstdev_b,2)],'Run b: -1std','Run b: Fit from Kp/Kd Circuits','Location','northwest')
legend('0.8nA input','0.8nA input+1std','0.8nA input-1std','Location','northwest')
% ylim([0 3])
xlim([0.5 100])
set(gcf, 'Color', 'w');
%% Base WN standard deviation
figure
semilogx(f_fftWNOL_a,magSys_a,'o')
hold on
semilogx(f_fftWNOL_a,magWN_a,'Linewidth',1.4,'Color','k')

semilogx(f_fftWNOL_a,magWN_a+PDstdev_a,'-.','Color','k')
semilogx(f_fftWNOL_a,magWN_a-PDstdev_a,'-.','Color','k')
% semilogx(f_fftWNOL_a,magKpKdWN_a,'--','Linewidth',1.4,'Color','k')
% semilogx(f_fftWNOL_b,magWN_b,'Linewidth',1.4,'Color','g')
% semilogx(f_fftWNOL_b,magWN_b+PDstdev_b,'-.','Color','g')
% semilogx(f_fftWNOL_b,magWN_b-PDstdev_b,'-.','Color','g')
% semilogx(f_fftWNOL_b,magKpKdWN_b,'--','Linewidth',1.4,'Color','g')
% semilogx(freq_Hz,magRatioOL_a,'-s','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_a+3*OLstd_a,'--','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_a-3*OLstd_a,'--','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_b,'-+','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_b+3*OLstd_b,'--','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_b-3*OLstd_b,'--','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_c,'-^','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_c+3*OLstd_c,'--','Linewidth',1.4)
% semilogx(freq_Hz,magRatioOL_c-3*OLstd_c,'--','Linewidth',1.4)
% fill(freq_Hz_flip,magflip,'b','FaceAlpha',0.1)
xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('White Noise Method Open Loop + 1 standard deviation')
legend('White Noise Data',['Run a: Fit Kp = ' num2str(KpPDFitWN_a,2) ',  Kd = ' num2str(KdPDFitWN_a,2)],['Run a: +1std = ' num2str(PDstdev_a,2)],'Run a: -1std','Location','northwest')
% legend(['Run a: Fit Kp = ' num2str(KpPDFitWN_a,2) ',  Kd = ' num2str(KdPDFitWN_a,2)],['Run b: Fit Kp = ' num2str(KpPDFitWN_b,2) ',  Kd = ' num2str(KdPDFitWN_b,2)],'Run a: Fit from Kp/Kd Circuits','Run b: Fit from Kp/Kd Circuits','0.1nA input','0.1nA input+3std','0.1nA input-3std','0.8nA input','0.8nA input+3std','0.8nA input-3std','4nA input','4nA input+3std','4nA input-3std','Location','northwest')
% ylim([0 4])
xlim([0.1 100])
set(gcf, 'Color', 'w');