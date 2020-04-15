% Design PD Controller for Inverted Pendulum
clear
clc

%Set up Peterka's experimental system
data = load('BVL_SurfaceTilt_1deg_EyesClosed.mat');
data = cell2struct(struct2cell(data), {'f','gain','phase','tranf'});   %rename tf within struct
f = [data(:).f];    %pull variables out of "data" struct
gain = [data(:).gain];
phase = [data(:).phase];
tranf = [data(:).tranf];
omega = 2*pi*f;
sys_exp = frd(tranf,omega);

% Peterka's "normal subject" constants
m = 83.3; %kg 
h = .896; %m - center of mass height
I = 81.1; %kg m^2
g = 9.81;% m/s

s = tf('s');

% Set up Inverted Pendulum TF
Gp = 1/(I*s^2-m*g*h);

gainLSR = 1; % amplifier gain applied to control signal (torque)
% gainFit = 5000/6;
gainFit = 833;

% parameters determined by eyeball fit
tauFit = 0.15;
KdFit = 0.05*6;   % for gain = 1, Kd = 250
KpFit = 0.21*6;   % for gain = 1, Kp = 1050

% parameters determined by least squares regression fit
KpLSR = 1007;   % Kp must be > mgh (732) for stability
KdLSR = 293;
tauLSR = 0.188;

% [tau,Kp,Kd] = LSRfitCLIP_MagPhase(gain2,f,gain,phase)

fb1 = 5;         % cutoff frequency in Hz (base 2)
lowFiltDerivative = (1/(s/(2*pi*fb1)+1));  % low pass filter

fb2 = 80;         % cutoff frequency in Hz (base 25)
lowFiltNeuron = (1/(s/(2*pi*fb2)+1));  % low pass filter

% Experimental Data
[expMag,expPhase] = bode(sys_exp,omega);
expMag = squeeze(expMag);
expMagdB = 20*log10(expMag);
expPhase = squeeze(expPhase);

% Hand Fit Model
[numFit,denFit] = pade(tauFit,1);    % pade delay approximation first order (in order to use impulse command)
delayFit = tf(numFit,denFit);
% delay = (-s+10)/(s+10);  

GcFit = KpFit + KdFit*s;
oltfFit = GcFit*Gp*gainFit;  % plus some delay

cltfFit = feedback(oltfFit,delayFit);  % G/(1+GH)
[fitMag,fitPhase] = bode(cltfFit,omega);
fitMag = squeeze(fitMag);
fitMagdB = 20*log10(fitMag);
fitPhase = squeeze(fitPhase);

[magCont,phaseCont,omegaCont] = bode(cltfFit);
magCont = squeeze(magCont);
magContdB = 20*log10(magCont);
phaseCont = squeeze(phaseCont);
fCont = omegaCont/(2*pi);

% LSR Fit Model
[numLSR,denLSR] = pade(tauLSR,1);    % pade delay approximation first order (in order to use impulse command)
delayLSR = tf(numLSR,denLSR);
% delay = (-s+10)/(s+10);  

GcLSR = KpLSR + KdLSR*s;
oltfLSR = GcLSR*Gp*gainLSR;  % plus some delay

cltfLSR = feedback(oltfLSR,delayLSR);  % G/(1+GH)
[LSRMag,LSRPhase] = bode(cltfLSR,omega);
LSRMag = squeeze(LSRMag);
LSRMagdB = 20*log10(LSRMag);
LSRPhase = squeeze(LSRPhase);

% no filter base PD fit model
data1 = load('tf2p1z.mat', 'tf2p1z');
cltf1 = [data1(:).tf2p1z];    %pull variables out of "data" struct
% cltf1 = feedback(tf1,1);
[cltfMag1,cltfPhase1] = bode(cltf1,omega);
cltfMag1 = squeeze(cltfMag1);
cltfPhase1 = squeeze(cltfPhase1);

% no filter base PD fit model no delay
data2 = load('tf2p1z_nodelay.mat','tf2p1z_nodelay');
cltf2 = [data2(:).tf2p1z_nodelay];    %pull variables out of "data" struct
% cltf2 = feedback(tf2,1);
[cltfMag2,cltfPhase2] = bode(cltf2,omega);
cltfMag2 = squeeze(cltfMag2);
cltfPhase2 = squeeze(cltfPhase2) - 360;

% % 2 80hz filter model
% data3 = load('CLTF_tf3_2filters_D_336.mat','tf3');
% tf3 = [data3(:).tf3];    %pull variables out of "data" struct
% cltf3 = feedback(tf3,1);
% [cltfMag3,cltfPhase3] = bode(cltf3,omega);
% cltfMag3 = squeeze(cltfMag3);
% cltfPhase3 = squeeze(cltfPhase3);
% 
% % 2 80hz filter model
% data4 = load('CLTF_tf4_1filters_D_224.mat','tf4');
% tf4 = [data4(:).tf4];    %pull variables out of "data" struct
% cltf4 = feedback(tf4,1);
% [cltfMag4,cltfPhase4] = bode(cltf4,omega);
% cltfMag4 = squeeze(cltfMag4);
% cltfPhase4 = squeeze(cltfPhase4);
%% Absolute Plots
figure
subplot(2,1,1)
semilogx(f,expMag,'-o')
hold on
% semilogx(f,fitMag)
% semilogx(f,LSRMag)
semilogx(fCont,magCont)
% semilogx(f,cltfMag1)
% semilogx(f,cltfMag2)
% semilogx(f,cltfMag3)
% semilogx(f,cltfMag4)
% legend('human data','4 filters','3 filters','2 filters','1 filter')
% legend('human data','my fit','Location','northeast')
title('Closed Loop Bode: Magnitude')
xlabel('Frequency Hz')
ylabel('Magnitude (abs)')
xlim([.01 10])
set(gcf, 'Color', 'w');

subplot(2,1,2)
semilogx(f,expPhase,'-o')
hold on
% semilogx(f,fitPhase)
% semilogx(f,LSRPhase)
semilogx(fCont,phaseCont)
% semilogx(f,cltfPhase1)
% semilogx(f,cltfPhase2)
% semilogx(f,cltfPhase3)
% semilogx(f,cltfPhase4)
% legend('human data','4 filters','3 filters','2 filters','1 filter')
legend('Human Data','Best Fit','Location','southwest')
title('Closed Loop Bode: Phase')
xlabel('Frequency Hz')
ylabel('Degrees')
xlim([.01 10])
set(gcf, 'Color', 'w');

%% dB Plots
figure
subplot(2,1,1)
semilogx(f,expMagdB,'-o')
hold on
semilogx(f,fitMagdB)
semilogx(f,LSRMagdB)
semilogx(fCont,magContdB)
title('Closed Loop Bode: Magnitude')
xlabel('Frequency Hz')
ylabel('Magnitude (dB)')
xlim([.01 10])
set(gcf, 'Color', 'w');

subplot(2,1,2)
semilogx(f,expPhase,'-o')
hold on
semilogx(f,fitPhase)
semilogx(f,LSRPhase)
semilogx(fCont,phaseCont)
% semilogx(f,cltfPhase1)
% semilogx(f,cltfPhase2)
% semilogx(f,cltfPhase3)
% semilogx(f,cltfPhase4)
% legend('human data','4 filters','3 filters','2 filters','1 filter')
legend('human data','my fit','LSR fit','Location','southwest')
title('Closed Loop Bode: Phase')
xlabel('Frequency Hz')
ylabel('Degrees')
xlim([.01 10])
set(gcf, 'Color', 'w');
% y = 2.86*impulse(cltf);    % 0.05 rads converted to nA
% 
% figure
% plot(y)
% ylabel('nA')
% 
% z = -1*impulse(cltf);     % 0.05 rad impulse (~2.9 degrees)
% 
% figure
% plot(z)
% ylabel('rad')