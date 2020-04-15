function [KpFitValue,stdev] = LSRfitKp_MagPhase_1Param_Error(Kpstart,Kpend,divs,freq_Hz,magRatioData,phaseData)
%Adapted from Cody Scharzenberger PSU's "LeastSquaresFittingScript.m" and
%included file AbsInd2DimInd.m
%6/5/2018

%Define the system parameter space domain in which to search.
% tau = linspace(0, 0.01, 10);
% tau = [0.004 1000];
Kp = linspace(Kpstart, Kpend, divs);
% Kp = [1.258 100];
% Kd = linspace(.1, 0.5, 15);
% Kd = [0.3 100];

s = tf('s');

fb1 = 5;         % cutoff frequency in Hz 
lowFiltDerivative = (1/(s/(2*pi*fb1)+1));  % low pass filter

fb2 = 80;         % cutoff frequency in Hz
lowFiltNeuron = (1/(s/(2*pi*fb2)+1));  % low pass filter

%Preallocate a variable to store the error values.
errors = zeros(length(Kp),1);

%Define the error weighting.
mag_weight = 1; phase_weight = 0.00;
% mag_weight = 0.95; phase_weight = 0.05;

%Compute the error associated with each parameter combination.
% for k1 = 1:length(tau)
    for k2 = 1:length(Kp)
%         for k3 = 1:length(Kd)
            
%             %Define the CL transfer function for PDIP system
%             [num,den] = pade(tau,1);    % pade delay approximation first order
%             delay = tf(num,den);
                   
%             oltf = delay*lowFiltNeuron^4*(Kp(k2) + Kd(k3)*s*lowFiltDerivative);
            oltf = lowFiltNeuron^2*Kp(k2);

            %Generate the frequency response data assocaited with this CL transfer function.
            [mag_est, phase_est] = bode(oltf, 2*pi*freq_Hz);
            
            %Reshape the frequency response data assocaited with this CL transfer function.
            [mag_est, phase_est] = deal( reshape(mag_est, size(freq_Hz)), reshape(phase_est, size(freq_Hz)) );
            
            %Compute the magnitude error.
            mag_error = norm(abs(mag_est - magRatioData));
            
            %Compute the phase error.
            phase_error = norm(abs(phase_est - phaseData));

            %Compute the weighted error.
            errors(k2) = mag_weight*mag_error + phase_weight*phase_error;
            
%         end
    end
% end

%Compute the minimum weighted error.
min_error = min(errors);

%Find the indexes associated with the minimum weighted error.
% inds = AbsInd2DimInd( size(errors), find(errors == min_error) );
inds = find(errors == min_error);

%Retrieve the natural frequency and damping ratio associated with minimum error.
KpFitValue = Kp(inds);

%Calculate standard deviation using fit value
oltfFit = lowFiltNeuron^2*KpFitValue;

[magFit, phaseFit] = bode(oltfFit, 2*pi*freq_Hz);
magFit = squeeze(magFit);
magFit = magFit';
magRatioData = magRatioData';

N = length(magRatioData);
var = magRatioData - magFit';
stdev = sqrt(sum(var.^2)/(N-1));

end


