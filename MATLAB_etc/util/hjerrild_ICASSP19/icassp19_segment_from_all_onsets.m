% ----------------------------------------------------------------------
% Onset detection and segmentation. 
%
%   INPUTS:
%           x:            observed recording
%           fs:           Sample rate to use
%           durationSec : Duration of each output segment (in seconds)
%   OUTPUTS:
%           segments:     a matrix containing all segments of the given
%                         duration extracted from the detected onset location.
%                         there is one segment for each onset.
%    onsetsInSeconds:     the onsets in seconds
% ------------------------------------------------------------------------------------------------
% segments = icassp19_segment_from_all_onsets(x,fs,durationSec)
% ------------------------------------------------------------------------------------------------
function [segments,onsetsInSeconds] = icassp19_segment_from_all_onsets(x,fs,durationSec)
x=[zeros(1000,1); x];
a = miraudio(x);
o=mironsets(a,'filter','diff','contrast',0.1);
seg = mirsegment(a,o);

onsetsInSeconds = mirgetdata(o);
onsetsInSamples = floor(onsetsInSeconds*fs)+1;

%% segment signal in to 40 ms durations
for n = 1: length(mirgetdata(o))
    segments(:,n) = x(onsetsInSamples(n):onsetsInSamples(n)+floor(durationSec*fs)+1);
    segments(:,n) = segments(:,n)/max(segments(:,n));
end
end

