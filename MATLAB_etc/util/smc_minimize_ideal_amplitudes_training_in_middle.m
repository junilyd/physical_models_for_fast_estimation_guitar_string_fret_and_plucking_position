% heuristic experiment that works well
%
function [pluckCmFromBridge, argmin, J, pCm,Cn] = smc_minimize_ideal_amplitudes_training_in_middle(amplitudes,L,trainingInMiddle,fn)
	if nargin < 4
		fn=1;
	end
	l = 1:length(amplitudes); % amplitude index
	amplitudes=amplitudes(l);
	
    if nargin<3, trainingInMiddle=0; end
    
	l(1)=[]; amplitudes(1)=[]; 
    h = 0.1; % pluck amplitude 
	Lhalf = L/2;
	pCm = 6:0.001:ceil(L/2); % search grid in cm.
	cnt=0;

    for ii=pCm % cm 
    	cnt=cnt+1;
    	R = ii/L;

    	Cn = ( 2*h )./( l.^2.*pi^2*R*(1-R) ) .* sin(l.*pi*R);
    	
        middleR = 1/2;
        if trainingInMiddle, Cn = l.*h.*abs(sin(l.*pi*ii/L))-( 2*h )./( l.^2.*pi^2*middleR*(1-middleR) ) .* sin(l.*pi*middleR), end
			
		J(cnt) = sum( (abs(diff(log(abs(Cn)))-diff(log(abs(amplitudes)))).^2));

	end

	[m,argmin] = min(J);
	pluckCmFromBridge = pCm(argmin);

	R=pCm(argmin)/L;
   	Cn = ( 2*h )./( l.^2.*pi^2*R*(1-R) ) .* sin(l.*pi*R);
end