% ----------------------------------------------------------------------
% Estimation of inharmonicity and pithc with inharmonic summation (approx. NLS)
%
%   INPUTS:
%           amplitudes  : Absolute value of partial amplitudes 
%           L           : Vibrating String Length (for placing P relative to the bridge) 
%   OUTPUTS:
%   pluckCmFromBridge   : Estimated plucking position PHat (in cm from bridge)
% ------------------------------------------------------------------------------------------------
% [pluckCmFromBridge] = icasssp19_plucking_position_estimator_LSD(amplitudes,L)
% ------------------------------------------------------------------------------------------------
function [pluckCmFromBridge] = icasssp19_plucking_position_estimator_LSD(amplitudes,L)
delta=1;
M = length(amplitudes);
m = (1:M)';
ii=0;
D_LS=[];
pCm = 6:0.001:ceil(L/2); % search grid in cm from the bridge.
for cnt=pCm
    ii=ii+1;
    P=cnt/L;
    Cm = 2*delta./(m.^2.*pi^2*P*(1-P)).*abs(sin((m.*pi*P))); 
    D_LS(ii) = sqrt( 1/M * sum( ( 10*log10(amplitudes(:)./Cm) ).^2) ) ;   
end
[~,argmin] = min(D_LS);
pluckCmFromBridge = pCm(argmin);
end