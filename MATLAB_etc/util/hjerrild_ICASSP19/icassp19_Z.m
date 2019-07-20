% ----------------------------------------------------------------
% Calculates Z-matrix containing the complex sinusoidal matrix representation. 
% In the proposed method it includes inharmonicity (piano model).
%
%    INPUT:
%           f0      Fundamental frequency in Hz. 
%           N       Given signal length
%           fs      Sample frequency
%           M       Number of partials in the model
%           B       Inharmonicity coefficient 
%                   (can be estimated using inharmonic summation)
%
%
%    OUTPUT
%           Z       Fourier matrix
%
% ----------------------------------
% function [Z] = icassp19_Z(f0,N,fs,M,B)
% ----------------------------------
%
function [Z] = icassp19_Z(f0,N,fs,M,B)
Z = zeros(N,M);
    for m = 1:M
        for n = 1:N
            Z(n,m) = exp((1i*2*pi*f0*m*sqrt(1+B*m^2))/fs*(n-1));
        end
    end
end