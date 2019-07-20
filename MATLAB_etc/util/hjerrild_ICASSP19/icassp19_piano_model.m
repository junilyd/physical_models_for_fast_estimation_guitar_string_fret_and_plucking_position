% -----------------------------------------------------------
% Calculate a harmonic frequency vector inharmonicity applied.
% (see Multi Pitch Estimation page 33)
%
%   Input 
%         beta  : stiffness coefficient (can be estimated from harmonic summation tuner) 
%         L     : number of harmonics that is desired
%   output
%         phi   : output vector with inharmonic frequencies
% -----------------------------------------
% function phi = smc_beta_model(f0,beta, L)
% -----------------------------------------
function phi = icassp19_piano_model(f0,beta, L)
    phi(1) = f0;
    l=2:L;
    phi(l) = f0*l.*sqrt(1+beta*l.^2);
end
