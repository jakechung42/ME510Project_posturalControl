function [f_fft,magRatioMod,phaseUnwrap,Ampin,Ampout] = bodebyFFT(inputSignal,outputSignal,sampleFrequency,numSegments,movAvgWinTime)
% Calculates absolute magnitude ratio and phase difference between input
% and output signals by fft(out)/fft(in)

dt = 1/sampleFrequency;
N = length(inputSignal);
pointsPerSegment = N/numSegments;
N1 = N/numSegments;
m1 = N1/2+1;
k = movAvgWinTime/dt;   % points for moving average window (time in s))

% Preprocess input/output signals to remove obvious large data errors and NaNs
x = inputSignal;
    indexNaNx = find(isnan(x));
    x(indexNaNx) = x(indexNaNx-1);
    indexBigx = find(x > 50);
    if indexBigx ~= 0
        fprintf('Outlier replaced in input signal: %d \n',x(indexBigx))
    end
    x(indexBigx) = 0;
y = outputSignal;
    indexNaNy = find(isnan(y));
    y(indexNaNy) = y(indexNaNy-1);
    indexBigy = find(y > 50);
    if indexBigy ~= 0
        fprintf('Outlier replaced in output signal: %d \n',y(indexBigy))
    end
    y(indexBigy) = 0;

xMat = reshape(x,[pointsPerSegment, numSegments])';
yMat = reshape(y,[pointsPerSegment, numSegments])';

% Note, the conversion factors are here for future reference for absolute
% calculations, but since they divide out in both the magnitude and phase
% calculation, the N/2 factor is unnecessary
FTx = fft(xMat,N1,2)/N1;
FTx = mean(FTx);
FTx = FTx(1:m1);    %single-sided spectrum
FTxPhase = FTx;
xthreshold = max(abs(FTx))/10000;
FTxPhase(abs(FTx) < xthreshold) = 0;

Ampin = abs(FTx); 
Ampin(2:end-1) = 2*Ampin(2:end-1);    %scale all but first and last value (amplitudes)
f_fft = linspace(0,sampleFrequency/2,m1);

FTy = fft(yMat,N1,2)/N1;
FTy = mean(FTy);
FTy = FTy(1:m1);    %single-sided spectrum
FTyPhase = FTy;
ythreshold = max(abs(FTy))/10000;
FTyPhase(abs(FTy) < ythreshold) = 0;

Ampout = abs(FTy);
Ampout(2:end-1) = 2*Ampout(2:end-1);    %scale all but first and last value

% Calculate Magnitude Ratio from fft(out)/fft(in) and Phase Diff between
% Input & Output Signals
magRat = Ampout./Ampin;
compRat = FTyPhase./FTxPhase;
phaseDiff = angle(compRat)*180/pi;  % wraps to 2pi by default - using unwrap makes data worse, can't really  use this phase to fit
% phaseDiff = atan2(imag(compRat),real(compRat))*180/pi; (same as angle
% function)

% unwrap phase (at least over first 10 Hz range)
phaseUnwrap = phaseDiff;
Mp = movmean(phaseUnwrap,k);
for jj = k/2+1:length(phaseUnwrap)
    diff = phaseUnwrap(jj) - Mp(jj-k/2);
    if diff >= 540
        phaseUnwrap(jj) = phaseUnwrap(jj) - 720;
    elseif diff >= 180
        phaseUnwrap(jj) = phaseUnwrap(jj) - 360;
    end
    Mp = movmean(phaseUnwrap,k);

end

%     figure
%     semilogx(f_fft,phaseDiff,'o')
%     hold on
%     semilogx(f_fft,phaseUnwrap,'x')
%     semilogx(f_fft,Mp)
%     xlim([0.01 10])

% remove magnitude outliers
Mm = movmean(magRat,k,'Endpoints','shrink');
Ms = movstd(magRat,k);

magRatioMod = magRat;

for ii = 2:length(magRatioMod)
    if magRatioMod(ii) >= Mm(ii) + Ms(ii)*2
        magRatioMod(ii) = magRatioMod(ii-1);
    elseif magRatioMod(ii) <= Mm(ii) - Ms(ii)*2
        magRatioMod(ii) = magRatioMod(ii-1);
    end
end

end

