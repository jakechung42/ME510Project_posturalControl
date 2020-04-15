% Compare bodes fitted to White Noise data for closed loop Controller overall, Kp and Kd circuits

clear
clc

% Import Data
% Closed Loop White Noise Data
% Run 1
D1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\NNinput_WNCL_70(1.5).dat');   % Desired input to NN
timeWNCL_a = D1(:,1);      % s
NNinput_desWNCL_a = D1(:,2);  

E1 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\NNangleOut_WNCL_70(1.5).dat');   
SysOutputWNCL_a = E1(:,2);  

A2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\ErrorScope_WNCL_70(1.5).dat');   
ErrorScopeWNCL_a = A2(:,2);  

B2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\Kp_WNCL_70(1.5).dat');   
KpScopeWNCL_a = B2(:,2); 

C2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\Kd_WNCL_70(1.5).dat');   
KdScopeWNCL_a = C2(:,2); 

D2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\kdpregain_WNCL_70(1.5).dat');   
kdpregainWNCL_a = D2(:,2); 

E2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\Kpscaler_WNCL_70(1.5).dat');   
KpscalerWNCL_a = E2(:,2); 

F2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\Kdscaler_WNCL_70(1.5).dat');   
KdscalerWNCL_a = F2(:,2); 

G2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\NNOutput_WNCL_70(1.5).dat');   
NNOutputWNCL_a = G2(:,2); 

H2 = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\ErrorN3_WNCL_70(1.5).dat');   
ErrorN3WNCL_a = H2(:,2); 

% Run 2
A = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\NNinput_WNCL_68(1.6).dat');   % Desired input to NN
timeWNCL_b = A(:,1);      % s
NNinput_desWNCL_b = A(:,2);  

B = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\NNangleOut_WNCL_68(1.6).dat');   
SysOutputWNCL_b = B(:,2);  

C = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\ErrorScope_WNCL_68(1.6).dat');   
ErrorScopeWNCL_b = C(:,2);  

D = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\Kp_WNCL_68(1.6).dat');   
KpScopeWNCL_b = D(:,2); 

E = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\Kd_WNCL_68(1.6).dat');   
KdScopeWNCL_b = E(:,2); 

F = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\kdpregain_WNCL_68(1.6).dat');   
kdpregainWNCL_b = F(:,2); 

G = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\Kpscaler_WNCL_68(1.6).dat');   
KpscalerWNCL_b = G(:,2); 

H = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\Kdscaler_WNCL_68(1.6).dat');   
KdscalerWNCL_b = H(:,2); 

I = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\NNOutput_WNCL_68(1.6).dat');   
NNOutputWNCL_b = I(:,2); 

J = importdata('C:\Users\Tiffany\Desktop\Sim_Outputs\WNCL_NNOutputs\ErrorN3_WNCL_68(1.6).dat');   
ErrorN3WNCL_b = J(:,2); 

%% Parameters

gainHuman = 833;
gainAB = 1388;
gainWN_a = 833;
gainWN_b = 1388;
relativeGainAdjust = gainWN_b/gainWN_a; % multiply b gains by this factor to compare to a for graphing
% gainWN = 833;
tau = 0.1;
tauHuman = 0.1;

s = tf('s');

% Peterka's "normal subject" constants
m = 83.3; %kg 
h = .896; %m - center of mass height
I = 81.1; %kg m^2
g = 9.81;% m/s

I_factor = 0.243734;    % assumes form of Inertia calculation is m*h^I_factor
m_adjust = 0;           % percent change to mass
h_adjust = 0;           % percent change to height

m2 = m + m_adjust*m;    %kg 
h2 = h + h_adjust*h;    %m - center of mass height
I2 = m*h^I_factor;      %kg m^2

Gp1 = 1/(I*s^2-m*g*h);
Gp2 = 1/(I2*s^2-m2*g*h2);

[num,den] = pade(tau,1);    % pade delay approximation first order
delay = tf(num,den);

[numHuman,denHuman] = pade(tauHuman,1);
delayHuman = tf(numHuman,denHuman);

%% Process White Noise Data - Run a

Gp = Gp1;
m = m;
h = h;
I = I;

dtWNCL_a = timeWNCL_a(3,1)-timeWNCL_a(2,1);
TsWNCL_a = 1/dtWNCL_a;
numSegmentsWNCL_a = 4;          % # segments to divide signal into for fft averaging
movAvgWinTimeWNCL_a = 0.5;      % window size for moving average in seconds

[f_fftWNCL_a,magSys_a,phaseSys_a,AmpinSys_a,AmpoutSys_a] = bodebyFFT(NNinput_desWNCL_a,SysOutputWNCL_a,TsWNCL_a,numSegmentsWNCL_a,movAvgWinTimeWNCL_a);
[f_fftWNCL_a,magKpgain_a,phaseKpgain_a,AmpinKpgain_a,AmpoutKpgain_a] = bodebyFFT(ErrorScopeWNCL_a,KpScopeWNCL_a,TsWNCL_a,numSegmentsWNCL_a,movAvgWinTimeWNCL_a);
[f_fftWNCL_a,magKdgain_a,phaseKdgain_a,AmpinKdgain_a,AmpoutKdgain_a] = bodebyFFT(kdpregainWNCL_a,KdScopeWNCL_a,TsWNCL_a,numSegmentsWNCL_a,movAvgWinTimeWNCL_a);

% Find avg measured Kp, Kd and standard deviations
mWindow_freq = 25;
mWindow_time = 5000;
kdbase = 0.046;  % measured from open loop
Kpmean_a = movmean(magKpgain_a(1:floor(end/50)),mWindow_freq);
Kpavgmean_a = mean(Kpmean_a)
Kpstd_a = movstd(magKpgain_a(1:floor(end/50)),mWindow_freq);
Kpavgstd_a = mean(Kpstd_a)
KpstdStatic_a = std(magKpgain_a(1:floor(end/50)))
Kdmean_a = movmean(magKdgain_a(1:floor(end/50)),mWindow_freq);
Kdavgmean_a = mean(Kdmean_a)*kdbase
Kdstd_a = movstd(magKdgain_a(1:floor(end/50)),mWindow_freq);
Kdavgstd_a = mean(Kdstd_a)*kdbase
KdstdStatic_a = std((magKdgain_a(1:floor(end/50)))*kdbase)
magKdgainstd_a = movstd(magSys_a(1:floor(end/50)),mWindow_freq);
KpscalerMean_a = movmean(KpscalerWNCL_a,mWindow_time);
KdscalerMean_a = movmean(KdscalerWNCL_a,mWindow_time);
ErrorMean_a = movmean(ErrorScopeWNCL_a,mWindow_time);
KpScopeMean_a = movmean(KpScopeWNCL_a,mWindow_time);
NNOutputMean_a = movmean(abs(NNOutputWNCL_a),mWindow_time);
ErrorN3Mean_a = movmean(abs(ErrorN3WNCL_a),mWindow_time);

figure
subplot(2,1,1)
semilogx(f_fftWNCL_a(1:floor(end/50)),Kpmean_a,'Linewidth',1.4)
hold on
semilogx(f_fftWNCL_a(1:floor(end/50)),Kpmean_a + Kpstd_a,'--','Linewidth',1.4)
semilogx(f_fftWNCL_a(1:floor(end/50)),max(Kpmean_a - Kpstd_a,0),'--','Linewidth',1.4)
semilogx(f_fftWNCL_a(1:floor(end/50)),magKpgain_a(1:floor(end/50)),'o')
title('Run a - Kp and kd gains by Frequency')
legend('Kp mean','Kp + 1stdev','Kp - 1stdev','data')
xlim([0.01 10])
set(gcf, 'Color', 'w');

subplot(2,1,2)
semilogx(f_fftWNCL_a(1:floor(end/50)),Kdmean_a*kdbase,'Linewidth',1.4)
hold on
semilogx(f_fftWNCL_a(1:floor(end/50)),Kdmean_a*kdbase + Kdstd_a*kdbase,'--','Linewidth',1.4)
semilogx(f_fftWNCL_a(1:floor(end/50)),max(Kdmean_a*kdbase - Kdstd_a*kdbase,0),'--','Linewidth',1.4)
semilogx(f_fftWNCL_a(1:floor(end/50)),magKdgain_a(1:floor(end/50))*kdbase,'o')
legend('Kd mean','Kd + 1stdev','Kd - 1stdev','data')
xlim([0.01 10])
set(gcf, 'Color', 'w');

% figure
% plot(timeWNCL_a,KpscalerWNCL_a,'o')
% hold on
% plot(timeWNCL_a,KpscalerMean_a,'Linewidth',1.4)
% title('Run a - Scaler')

% figure
% plot(timeWNCL_a,KdscalerWNCL_a,'o')
% hold on
% plot(timeWNCL_a,KdscalerMean_a,'Linewidth',1.4)
% title('Kd Scaler')

% figure
% plot(timeWNCL_a,ErrorScopeWNCL_a)
% hold on
% plot(timeWNCL_a,ErrorMean_a,'Linewidth',2)
% plot(timeWNCL_a,KpScopeWNCL_a)
% plot(timeWNCL_a,KpScopeMean_a,'Linewidth',2)
% title('Kp Gain over Time')
% legend('Error Scope','Error Mean','Kp Scope','Kp Mean')

% figure
% plot(timeWNCL_a,ErrorN3WNCL_a)
% hold on
% plot(timeWNCL_a,ErrorN3Mean_a,'Linewidth',2)
% plot(timeWNCL_a,NNOutputWNCL_a)
% plot(timeWNCL_a,NNOutputMean_a,'Linewidth',2)
% title('Error and Torque over Time')
% legend('Error','Error Mean','Output Torque','Torque Mean')

% Fit Closed Loop WN free fit (will probably be inaccurate)
Kpstart = 0.1; Kpend = 2.5; Kpdivs = 40;
Kdstart = 0.1; Kdend = 1; Kddivs = 40;
[KpFitWN_a,KdFitWN_a,LSRFit_a_stdev] = LSRfitCLPD_MagPhase_2Param_Error_VarySys(Kpstart,Kpend,Kpdivs,Kdstart,Kdend,Kddivs,f_fftWNCL_a(1:floor(end/50)),tau,gainWN_a,magSys_a(1:floor(end/50)),phaseSys_a(1:floor(end/50)),m,h,I)

GcWN_a = KpFitWN_a + KdFitWN_a*s;
oltfWN_a = delay*GcWN_a*Gp*gainWN_a;
cltfWN_a = feedback(oltfWN_a,1);

f_fftWNCL_rad_a = f_fftWNCL_a*2*pi;
[magWN_a,phaseWN_a] = bode(cltfWN_a,f_fftWNCL_rad_a);    % requires frequency in rad/s
magWN_a = squeeze(magWN_a);
phaseWN_a = squeeze(phaseWN_a);

% Closed loop TFs from Measured Kp/Kd and stdevs
% with mean Kp/Kd
GcMean_a = Kpavgmean_a + Kdavgmean_a*s;
oltfMean_a = delay*GcMean_a*Gp*gainWN_a;
cltfMean_a = feedback(oltfMean_a,1);

[magMean_a,phaseMean_a] = bode(cltfMean_a,f_fftWNCL_rad_a);    % requires frequency in rad/s
magMean_a = squeeze(magMean_a);
phaseMean_a = squeeze(phaseMean_a);

% with high Kp/Kd (mean + 1std)
KpHigh_a = Kpavgmean_a + Kpavgstd_a;
KdHigh_a = Kdavgmean_a + Kdavgstd_a;
GcHigh_a = KpHigh_a + KdHigh_a*s;
oltfHigh_a = delay*GcHigh_a*Gp*gainWN_a;
cltfHigh_a = feedback(oltfHigh_a,1);

[magHigh_a,phaseHigh_a] = bode(cltfHigh_a,f_fftWNCL_rad_a);    % requires frequency in rad/s
magHigh_a = squeeze(magHigh_a);
phaseHigh_a = squeeze(phaseHigh_a);

KpLow_a = max([(Kpavgmean_a - Kpavgstd_a) 0])
KdLow_a = max([(Kdavgmean_a - Kdavgstd_a) 0])
GcLow_a = KpLow_a + KdLow_a*s;
oltfLow_a = delay*GcLow_a*Gp*gainWN_a;
cltfLow_a = feedback(oltfLow_a,1);

[magLow_a,phaseLow_a] = bode(cltfLow_a,f_fftWNCL_rad_a);    % requires frequency in rad/s
magLow_a = squeeze(magLow_a);
phaseLow_a = squeeze(phaseLow_a);

%% Process White Noise Data - Run b

Gp = Gp2;
m = m2;
h = h2;
I = I2;

dtWNCL_b = timeWNCL_b(3,1)-timeWNCL_b(2,1);
TsWNCL_b = 1/dtWNCL_b;
numSegmentsWNCL_b = 4;   % # segments to divide signal into for fft averaging
movAvgWinTimeWNCL_b = 0.5;   % window size for moving average in seconds

[f_fftWNCL_b,magSys_b,phaseSys_b,AmpinSys_b,AmpoutSys_b] = bodebyFFT(NNinput_desWNCL_b,SysOutputWNCL_b,TsWNCL_b,numSegmentsWNCL_b,movAvgWinTimeWNCL_b);
[f_fftWNCL_b,magKpgain_b,phaseKpgain_b,AmpinKpgain_b,AmpoutKpgain_b] = bodebyFFT(ErrorScopeWNCL_b,KpScopeWNCL_b,TsWNCL_b,numSegmentsWNCL_b,movAvgWinTimeWNCL_b);
[f_fftWNCL_b,magKdgain_b,phaseKdgain_b,AmpinKdgain_b,AmpoutKdgain_b] = bodebyFFT(kdpregainWNCL_b,KdScopeWNCL_b,TsWNCL_b,numSegmentsWNCL_b,movAvgWinTimeWNCL_b);

% Find avg measured Kp, Kd and standard deviations
mWindow_freq = 25;
mWindow_time = 5000;
kdbase = 0.046;  % measured from open loop
Kpmean_b = movmean(magKpgain_b(1:floor(end/50)),mWindow_freq);
Kpavgmean_b = mean(Kpmean_b)
Kpavgmean_b_adjust = Kpavgmean_b*relativeGainAdjust
Kpstd_b = movstd(magKpgain_b(1:floor(end/50)),mWindow_freq);
Kpavgstd_b = mean(Kpstd_b);
Kpavgstd_b_adjust = Kpavgstd_b*relativeGainAdjust
KpstdStatic_b = std(magKpgain_b(1:floor(end/50)))
Kdmean_b = movmean(magKdgain_b(1:floor(end/50)),mWindow_freq);
Kdavgmean_b = mean(Kdmean_b)*kdbase
Kdavgmean_b_adjust = Kdavgmean_b*relativeGainAdjust
Kdstd_b = movstd(magKdgain_b(1:floor(end/50)),mWindow_freq);
Kdavgstd_b = mean(Kdstd_b)*kdbase
Kdavgstd_b_adjust = Kdavgstd_b*relativeGainAdjust
KdstdStatic_b = std((magKdgain_b(1:floor(end/50)))*kdbase)
magKdgainstd_b = movstd(magSys_b(1:floor(end/50)),mWindow_freq);
KpscalerMean_b = movmean(KpscalerWNCL_b,mWindow_time);
KdscalerMean_b = movmean(KdscalerWNCL_b,mWindow_time);
ErrorMean_b = movmean(ErrorScopeWNCL_b,mWindow_time);
KpScopeMean_b = movmean(KpScopeWNCL_b,mWindow_time);
NNOutputMean_b = movmean(abs(NNOutputWNCL_b),mWindow_time);
ErrorN3Mean_b = movmean(abs(ErrorN3WNCL_b),mWindow_time);

figure
subplot(2,1,1)
semilogx(f_fftWNCL_b(1:floor(end/50)),magKpgain_b(1:floor(end/50)),'o')
hold on
semilogx(f_fftWNCL_b(1:floor(end/50)),Kpmean_b,'Linewidth',1.4,'Color','k')
semilogx(f_fftWNCL_b(1:floor(end/50)),Kpmean_b + Kpstd_b,'--','Linewidth',1.4,'Color','k')
semilogx(f_fftWNCL_b(1:floor(end/50)),max(Kpmean_b - Kpstd_b,0),'--','Linewidth',1.4,'Color','k')
title('Run b - Kp and Kd gains')
legend('Kp data','Kp mean','Kp + 1stdev','Kp - 1stdev')
xlim([0.01 10])
ylim([0.73 0.77])
set(gcf, 'Color', 'w');
% xlabel('Frequency (Hz)')
ylabel('Kp Gain')

subplot(2,1,2)
semilogx(f_fftWNCL_b(1:floor(end/50)),magKdgain_b(1:floor(end/50))*kdbase,'o')
hold on
semilogx(f_fftWNCL_b(1:floor(end/50)),Kdmean_b*kdbase,'Linewidth',1.4,'Color','k')
semilogx(f_fftWNCL_b(1:floor(end/50)),Kdmean_b*kdbase + Kdstd_b*kdbase,'--','Linewidth',1.4,'Color','k')
semilogx(f_fftWNCL_b(1:floor(end/50)),max(Kdmean_b*kdbase - Kdstd_b*kdbase,0),'--','Linewidth',1.4,'Color','k')
legend('Kd data','Kd mean','Kd + 1stdev','Kd - 1stdev')
xlim([0.01 10])
ylim([0.185 0.20])
set(gcf, 'Color', 'w');
xlabel('Frequency (Hz)')
ylabel('Kd Gain')

figure
plot(timeWNCL_b,KpscalerWNCL_b,'o')
hold on
plot(timeWNCL_b,KpscalerMean_b,'Linewidth',1.4)
title('Run b - Scaler')

% figure
% plot(timeWNCL_b,KdscalerWNCL_b,'o')
% hold on
% plot(timeWNCL_b,KdscalerMean_b,'Linewidth',1.4)
% title('Kd Scaler')

% figure
% plot(timeWNCL_b,ErrorScopeWNCL_b)
% hold on
% plot(timeWNCL_b,ErrorMean_b,'Linewidth',2)
% plot(timeWNCL_b,KpScopeWNCL_b)
% plot(timeWNCL_b,KpScopeMean_b,'Linewidth',2)
% title('Run b - Kp Gain over Time')
% legend('Error Scope','Error Mean','Kp Scope','Kp Mean')
% 
% figure
% plot(timeWNCL_b,ErrorN3WNCL_b)
% hold on
% plot(timeWNCL_b,ErrorN3Mean_b,'Linewidth',2)
% plot(timeWNCL_b,NNOutputWNCL_b)
% plot(timeWNCL_b,NNOutputMean_b,'Linewidth',2)
% title('Run b - Error and Torque over Time')
% legend('Error','Error Mean','Output Torque','Torque Mean')

% Fit Closed Loop WN free (will probably be inaccurate)
Kpstart = 0.1; Kpend = 2.5; Kpdivs = 40;
Kdstart = 0.1; Kdend = 1; Kddivs = 40;
[KpFitWN_b,KdFitWN_b,LSRFit_b_stdev] = LSRfitCLPD_MagPhase_2Param_Error_VarySys(Kpstart,Kpend,Kpdivs,Kdstart,Kdend,Kddivs,f_fftWNCL_b(1:floor(end/50)),tau,gainWN_b,magSys_b(1:floor(end/50)),phaseSys_b(1:floor(end/50)),m,h,I);

GcWN_b = KpFitWN_b + KdFitWN_b*s;
oltfWN_b = delay*GcWN_b*Gp*gainWN_b;
cltfWN_b = feedback(oltfWN_b,1);

KpFitWN_b_adjust = KpFitWN_b*relativeGainAdjust
KdFitWN_b_adjust = KdFitWN_b*relativeGainAdjust

f_fftWNCL_rad_b = f_fftWNCL_b*2*pi;
[magWN_b,phaseWN_b] = bode(cltfWN_b,f_fftWNCL_rad_b);    % requires frequency in rad/s
magWN_b = squeeze(magWN_b);
phaseWN_b = squeeze(phaseWN_b);

% Closed loop TFs from Measured Kp/Kd and stdevs
% with mean Kp/Kd
GcMean_b = Kpavgmean_b + Kdavgmean_b*s;
oltfMean_b = delay*GcMean_b*Gp*gainWN_b;
cltfMean_b = feedback(oltfMean_b,1);

[magMean_b,phaseMean_b] = bode(cltfMean_b,f_fftWNCL_rad_b);    % requires frequency in rad/s
magMean_b = squeeze(magMean_b);
phaseMean_b = squeeze(phaseMean_b);

% with high Kp/Kd (mean + 1std)
KpHigh_b = Kpavgmean_b + Kpavgstd_b
KdHigh_b = Kdavgmean_b + Kdavgstd_b
GcHigh_b = KpHigh_b + KdHigh_b*s;
oltfHigh_b = delay*GcHigh_b*Gp*gainWN_b;
cltfHigh_b = feedback(oltfHigh_b,1);

[magHigh_b,phaseHigh_b] = bode(cltfHigh_b,f_fftWNCL_rad_b);    % requires frequency in rad/s
magHigh_b = squeeze(magHigh_b);
phaseHigh_b = squeeze(phaseHigh_b);

KpLow_b = max([(Kpavgmean_b - Kpavgstd_b) 0])
KdLow_b = max([(Kdavgmean_b - Kdavgstd_b) 0])
GcLow_b = KpLow_b + KdLow_b*s;
oltfLow_b = delay*GcLow_b*Gp*gainWN_b;
cltfLow_b = feedback(oltfLow_b,1);

[magLow_b,phaseLow_b] = bode(cltfLow_b,f_fftWNCL_rad_b);    % requires frequency in rad/s
magLow_b = squeeze(magLow_b);
phaseLow_b = squeeze(phaseLow_b);

%% Closed Loop Human Reference Model & Auto Bode Fits

% Fit TF to Human Data (STATIC)
Kd = 0.05*6;   % for gain = 1, Kd = 250
Kp = 0.21*6;   % for gain = 1, Kp = 1050

% fb1 = 5;         % cutoff frequency in Hz (base 2)
% lowFiltDerivative = (1/(s/(2*pi*fb1)+1));  % low pass filter
% 
% fb2 = 80;         % cutoff frequency in Hz (base 80)
% lowFiltNeuron = (1/(s/(2*pi*fb2)+1));  % low pass filter

Gc = Kp + Kd*s;
oltf = delayHuman*Gc*Gp*gainHuman;

cltf = feedback(oltf,1);

[magHuman,phaseHuman,omegaHuman] = bode(cltf);    % requires frequency in rad/s
magHuman = squeeze(magHuman);
phaseHuman = squeeze(phaseHuman);
f_Human = omegaHuman/(2*pi);

%% Plots

% CLosed Loop figure - Run a
figure
semilogx(f_Human,magHuman,'Linewidth',1.2)
hold on
semilogx(f_fftWNCL_a,magSys_a,'o')
semilogx(f_fftWNCL_a,magWN_a,'-.','Linewidth',1.4)
semilogx(f_fftWNCL_a,magMean_a,'Linewidth',1.4)
semilogx(f_fftWNCL_a,magHigh_a,'--','Linewidth',1.4)
semilogx(f_fftWNCL_a,magLow_a,'--','Linewidth',1.4)
xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('Closed Loop System Transfer Function with Variable Inputs (Run a)')
legend('Human Fit','White Noise Data','White Noise Free Fit','Measured Mean Kp/Kd Fit','High Fit (Kp/Kd + 1std)','Low Fit (Kp/Kd - 1std)','Location','northwest')
% ylim([-40 20])
xlim([0.01 10])
set(gcf, 'Color', 'w');
% export_fig('C:\Users\Tiffany\Desktop\Sim_Outputs\Figs\CLABWN_MultiInput_Weighting1_WN26', '-bmp','-q101','-nocrop','-a1');

% CLosed Loop figure - Run b
figure
semilogx(f_Human,magHuman,'Linewidth',1.2)
hold on
semilogx(f_fftWNCL_b,magSys_b,'o')
semilogx(f_fftWNCL_b,magWN_b,'-.','Linewidth',1.4)
semilogx(f_fftWNCL_b,magMean_b,'Linewidth',1.4)
semilogx(f_fftWNCL_b,magHigh_b,'--','Linewidth',1.4)
semilogx(f_fftWNCL_b,magLow_b,'--','Linewidth',1.4)
xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('Closed Loop System Transfer Function with Variable Inputs (Run b)')
legend('Human Fit','White Noise Data','White Noise Free Fit','Measured Mean Kp/Kd Fit','High Fit (Kp/Kd + 1std)','Low Fit (Kp/Kd - 1std)','Location','northeast')
% ylim([-40 20])
xlim([0.01 10])
set(gcf, 'Color', 'w');
% export_fig('C:\Users\Tiffany\Desktop\Sim_Outputs\Figs\CLABWN_MultiInput_Weighting1_WN26', '-bmp','-q101','-nocrop','-a1');

% CLosed Loop figure - Compare
figure
semilogx(f_Human,magHuman,'Linewidth',1.2)
hold on
semilogx(f_fftWNCL_a,magSys_a,'o')
semilogx(f_fftWNCL_a,magWN_a,'-.','Linewidth',1.4)
semilogx(f_fftWNCL_a,magMean_a,'Linewidth',1.4)
semilogx(f_fftWNCL_a,magHigh_a,'--','Linewidth',1.4)
semilogx(f_fftWNCL_a,magLow_a,'--','Linewidth',1.4)
semilogx(f_fftWNCL_b,magSys_b,'x')
semilogx(f_fftWNCL_b,magWN_b,'-.','Linewidth',1.4)
semilogx(f_fftWNCL_b,magMean_b,'Linewidth',1.4)
semilogx(f_fftWNCL_b,magHigh_b,'--','Linewidth',1.4)
semilogx(f_fftWNCL_b,magLow_b,'--','Linewidth',1.4)
xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('Closed Loop System Transfer Function with Variable Inputs Compare')
legend('Human Fit','White Noise Data Run a','White Noise Free Fit Run a','Measured Mean Kp/Kd Fit Run a','High Fit (Kp/Kd + 1std) Run a','Low Fit (Kp/Kd - 1std) Run a','White Noise Data Run b','White Noise Free Fit Run b','Measured Mean Kp/Kd Fit Run b','High Fit (Kp/Kd + 1std) Run b','Low Fit (Kp/Kd - 1std) Run b','Location','northeast')
% ylim([-40 20])
xlim([0.01 10])
set(gcf, 'Color', 'w');
% export_fig('C:\Users\Tiffany\Desktop\Sim_Outputs\Figs\CLABWN_Compare_Premod11a_Postmod32', '-bmp','-q101','-nocrop','-a1');

% CLosed Loop figure - Compare
figure
% semilogx(f_Human,magHuman,'Linewidth',1.2)
semilogx(f_fftWNCL_a,magSys_a,'o','Color','k')
hold on
semilogx(f_fftWNCL_a,magWN_a,'-','Linewidth',1.4,'Color','k')
% semilogx(f_fftWNCL_a,magMean_a,'-.','Linewidth',1.4)
% semilogx(f_fftWNCL_a(27:floor(end/50)),magWN_a_High(27:end),'--','Linewidth',1.4)
% semilogx(f_fftWNCL_a(27:floor(end/50)),magWN_a_Low(27:end),'--','Linewidth',1.4)
semilogx(f_fftWNCL_b,magSys_b,'x','Color','m')
semilogx(f_fftWNCL_b,magWN_b,'-','Linewidth',1.4,'Color','m')
% semilogx(f_fftWNCL_b,magMean_b,'-.','Linewidth',1.4)
% semilogx(f_fftWNCL_b(27:floor(end/50)),magWN_b_High(27:end),'--','Linewidth',1.4)
% semilogx(f_fftWNCL_b(27:floor(end/50)),magWN_b_Low(27:end),'--','Linewidth',1.4)
% fill(freq_Hz_flip,magflip,'b','FaceAlpha',0.1)
xlabel('Frequency (Hz)')
ylabel('Mag Ratio (abs)')
title('(b) Closed Loop System: Input 0.6nA, Sensitivity Test - Inertia 152%')
% legend('Run a: White Noise Data',['Run a: White Noise Fit, Kp = ' num2str(KpFitWN_a,2) ' Kd = ' num2str(KdFitWN_a,2) ' Std = ' num2str(LSRFit_a_stdev,2)],['Run a: Fit from Avg Measured Kp = ' num2str(Kpavgmean_a,2) ' +- ' num2str(KpstdStatic_a,2) ',  Kd = ' num2str(Kdavgmean_a,2) ' +- ' num2str(KdstdStatic_a,2) ],...
%         'Run b: White Noise Data ',['Run b: White Noise Fit, Kp = ' num2str(KpFitWN_b,2) ' Kd = ' num2str(KdFitWN_b,2) ' Std = ' num2str(LSRFit_b_stdev,2)],['Run b: Fit from Avg Measured Kp = ' num2str(Kpavgmean_b,2) ' +- ' num2str(KpstdStatic_b,2)  ',  Kd = ' num2str(Kdavgmean_b,2) ' +- ' num2str(KdstdStatic_b,2)],'Location','northeast')
legend('Baseline Inertia: White Noise Data',['Baseline Inertia: White Noise Fit, Kp = ' num2str(KpFitWN_a,2) ' Kd = ' num2str(KdFitWN_a,2) ' Std = ' num2str(LSRFit_a_stdev,2)],...
        'Inertia 152%: White Noise Data ',['Inertia 152%: White Noise Fit, Kp = ' num2str(KpFitWN_b,2) ' Kd = ' num2str(KdFitWN_b,2) ' Std = ' num2str(LSRFit_b_stdev,2)],'Location','northeast')
ylim([0 14])
xlim([0.01 3])
set(gcf, 'Color', 'w');
% export_fig('C:\Users\Tiffany\Desktop\Sim_Outputs\Figs\CLABWN_Compare_Postmod32_Weight34_datastd', '-bmp','-q101','-nocrop','-a1');

