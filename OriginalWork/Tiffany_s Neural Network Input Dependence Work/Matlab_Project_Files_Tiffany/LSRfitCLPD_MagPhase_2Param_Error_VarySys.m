function [KpFitValue,KdFitValue,stdev] = LSRfitCLPD_MagPhase_2Param_VarySys(Kpstart,Kpend,Kpdivs,Kdstart,Kdend,Kddivs,freq_Hz,tau,gain,magRatioData,phaseData,m,h,I)
%Adapted from Cody Scharzenberger PSU's "LeastSquaresFittingScript.m" and
%included file AbsInd2DimInd.m
%6/5/2018

%Define the system parameter space domain in which to search.
% tau = linspace(0, 0.01, 10);
% tau = [0.004 1000];
Kp = linspace(Kpstart, Kpend, Kpdivs);
Kd = linspace(Kdstart, Kdend, Kddivs);

% Peterka's "normal subject" constants
% m = 83.3; %kg 
% h = .896; %m - center of mass height
% I = 81.1; %kg m^2
g = 9.81;% m/s

s = tf('s');
Gp = 1/(I*s^2-m*g*h);

fb1 = 5;         % cutoff frequency in Hz 
lowFiltDerivative = (1/(s/(2*pi*fb1)+1));  % low pass filter

fb2 = 80;         % cutoff frequency in Hz
lowFiltNeuron = (1/(s/(2*pi*fb2)+1));  % low pass filter

%Preallocate a variable to store the error values.
errors = zeros(length(Kd), length(Kp));

%Define the error weighting.
mag_weight = 1; phase_weight = 0;
% mag_weight = 0.99; phase_weight = 0.01;

%Compute the error associated with each parameter combination.
% for k1 = 1:length(tau)
    for k2 = 1:length(Kp)
        for k3 = 1:length(Kd)
            
            %Define the CL transfer function for PDIP system
            [num,den] = pade(tau,1);    % pade delay approximation first order
            delay = tf(num,den);
                   
            Gc = Kp(k2) + Kd(k3)*s;
            oltf = Gc*Gp*gain*delay;
            cltf = feedback(oltf,1);

            %Generate the frequency response data assocaited with this CL transfer function.
            [mag_est, phase_est] = bode(cltf, 2*pi*freq_Hz);
            
            %Reshape the frequency response data assocaited with this CL transfer function.
            [mag_est, phase_est] = deal( reshape(mag_est, size(freq_Hz)), reshape(phase_est, size(freq_Hz)) );
            
            %Compute the magnitude error.
            mag_error = norm(abs(mag_est - magRatioData));
            
            %Compute the phase error.
            phase_error = norm(abs(phase_est - phaseData));

            %Compute the weighted error.
            errors(k3, k2) = mag_weight*mag_error + phase_weight*phase_error;
            
        end
    end
% end

%Compute the minimum weighted error.
min_error = min(min(errors));

%Find the indexes associated with the minimum weighted error.
inds = AbsInd2DimInd( size(errors), find(errors == min_error) );
% inds = find(errors == min_error);

%Retrieve the natural frequency and damping ratio associated with minimum error.
KdFitValue = Kd(inds(1)); KpFitValue = Kp(inds(2));

%Calculate standard deviation using fit value
Gc = KpFitValue + KdFitValue*s;
oltf = Gc*Gp*gain*delay;
cltf = feedback(oltf,1);

[magFit, phaseFit] = bode(cltf, 2*pi*freq_Hz);
magFit = squeeze(magFit);
magFit = magFit';
magRatioData = magRatioData';

N = length(magRatioData);
var = magRatioData - magFit';
stdev = sqrt(sum(var.^2)/(N-1));

end


