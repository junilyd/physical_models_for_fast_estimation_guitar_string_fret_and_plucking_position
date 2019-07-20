% ----------------------------------------------------------------------
% This function creates an index vector for signal partials 
% used for harmonic summation.
% The fundamental frequency input vector and the FFT-length, 
% will create the dimensions of the output.
% The peaks are placed in a zero vector, on the correct frequency axis.
%
%   INPUTS:     
%               fft_sig:    fft-signal which is real valued. This is assumed to be cut at  
%                           Nyquist frequencu. (1:fs/2)
%               fs:         sample frequency used for FFT.
%               L:          The desired number of harmonics
%               f0Area:     frequency searc area
%
%   OUTPUTS:
%               index:    Index of L inharmonics
%
% ---------------------------------------------------------------------
% function [index] = smc_inharmonic_index(fft_sig, fs, L, f0Area,beta)
% ---------------------------------------------------------------------
function [index] = icassp19_inharmonic_index(fft_sig, fs, L, f0Area, beta)
l=1:L;
N = length(fft_sig);

phi = (icassp19_piano_model(f0Area, beta, L ))';

index_ratio = N/(fs/2);
        index = round(phi(l).*(2*N/fs)+1);
        index=index(:);
end
