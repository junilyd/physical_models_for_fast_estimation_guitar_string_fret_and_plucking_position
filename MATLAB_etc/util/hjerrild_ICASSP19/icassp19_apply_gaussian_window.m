% -------------------------------------------------------------------
% Apply a Gaussian window.
% 
% INPUT: 
%         x               :  The signal which should is wanted to be
%                            modified
%         a               :  width/variance control. (optional)
%						  :  inversly proportional to the standard deviation of a normal Gaussian distro.		
% OUTPUT:
%         x	              :  windowed signal
%
% ---------------------------------------------------
% x = icassp19_apply_gaussian_window(x)
% ---------------------------------------------------
function x = icassp19_apply_gaussian_window(x,a)
       	Nwin = length(x)-1;
       	if nargin < 2, 
       		a = 2.5; 
        end
       	n = (0:Nwin)' - Nwin/2;
		w = exp(-(1/2)*(a*n/(Nwin/2)).^2);
		x = x.*w;
end