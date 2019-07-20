%   ----------------------------------------------------------
%   Estimates the pitch from harmonic summation
%   Here the tuner also zooms in to get some calc. efficiency
%   --------------------------------------------------------------------------------
%   INPUTS:
%           X:       fft of input signal
%           f0Limits: Vector of fundamental frequrncies in hertz (the search area)
%           M:       Number of Harmonics to use.
%           fs:      Sampling rate to use (usually samplerate of input signal)
%   OUTPUTS:
%           pitch:   The pitch estimate
% --------------------------------------------------------
% pitch3 = icassp19_harmonic_summation(X, f0Limits, M, fs)
% --------------------------------------------------------
function [pitch3,f0Index] = icassp19_harmonic_summation(X, f0Limits, M, fs)

f0Grid = [min(f0Limits):0.1:max(f0Limits)];

i=1;
for f=f0Grid
    [index] = icassp19_harmonic_index(X, fs, M, f);
    cost(i) = sum(X(index)); i = i+1;
end
[C, I] = max(cost);
pitch = f0Grid(I);

f0_area2 = [pitch-4:0.01:pitch+4];
i=1;
for f=f0_area2
    [index] = icassp19_harmonic_index(X, fs, M, f);
    cost2(i) = sum(X(index));i=i+1;
end
[C,I] = max(cost2);
pitch2 = f0_area2(I);

f0_area3 = [pitch2-0.01:0.001:pitch2+0.01];
i=1;
for f=f0_area3
    [index] = icassp19_harmonic_index(X, fs, M, f);
    cost3(i) = sum(X(index));i=i+1;
end
[C,I] = max(cost3);
pitch3 = f0_area3(I);

% return maximum near f0 estimate
[firstPeakIndex] = icassp19_harmonic_index(X, fs, 1, pitch3);
[lowerIndex] = icassp19_harmonic_index(X, fs, 1, pitch3-pitch3/4);
[upperIndex] = icassp19_harmonic_index(X, fs, 1, pitch3+pitch3/4);

[~,peakNdx] = max(X(lowerIndex:upperIndex));
f0Index=peakNdx+lowerIndex-1;


end
