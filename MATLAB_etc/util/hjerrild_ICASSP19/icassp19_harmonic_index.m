% ----------------------------------------------------------------------
%  This function creates a vector with indeces of the partials in the signal.
%   The indeces can then be used for sparse harmonic summation
% 
%   INPUTS:     
%               fft_sig:    fft-signal which is real valued. This is assumed to be cut at  
%                           Nyquist frequencu. (1:fs/2)
%               fs:         sample frequency used for FFT.
%               L:          The desired number of harmonics
%               f0Area:     Vector containg the search area
%
%   OUTPUTS:
%               index:    Index of L harmonics
%
% --------------------------------------------------------------
% function [index] = smc_harmonic_index(fft_sig, fs, L, f0Area)
% --------------------------------------------------------------
function [index] = icassp19_harmonic_index(fft_sig, fs, L, f0Area)
    l=1:L;
    f0Area = f0Area(:);
    NFFT = length(fft_sig);
    f = f0Area*l;
    index = round(f(l).*(2*NFFT/fs)+1);
    index = index(:);
end
