% ----------------------------------------------------------------------
% Calculate analytic complex signal from real signal   
%   Calcualates the analytic signal from a real signal into size of next
%   power of 2 from lenght(sig).
%
%   INPUT: 
%           sig: real signal
%   OUTPUT:
%           x: analytic signal (complex signal)
%
% --------------------
% x = icassp19_hilbert_transform(sig,fs)
% --------------------
function x = icassp19_hilbert_transform(sig, fs)

N = length(sig);
if mod(length(sig),2) == 1
   sig = sig(1:end-1);
   N = N-1;
end
fft_long = fft(sig, N);
X = 2*abs(fft_long(1:N/2)).^2; 
f_axis   = fs/2 * linspace(0,1,N/2);
h = zeros(1,N)';
h(1) = 1;
h(2:N/2) = 2;
h(N/2+1) = 1;
x_fft = zeros(1,N)';
for i = 1:N
    x_fft(i) = fft_long(i)*h(i);
end
x = ifft(x_fft);
end
