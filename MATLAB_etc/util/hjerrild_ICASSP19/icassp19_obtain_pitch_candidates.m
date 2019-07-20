%   --------------------------------------------------------
%   Select candidate all classes (string and fret combinations) based on the equal
%   tempered scale and the trained model of pitch for the given classes. 
%   --------------------------------------------------------
%   INPUTS:
%   observedPitch      :	Observed pitch estimate in Hz 
%	frequencyReference :    trained model of pitch for the given classes.
%
%   OUTPUTS:
%       strind and fret candidates for the given pitch.
%		(other string and fret combinations can have a 
%		probability of 0. (i.e., see abesser et al. [7])
% 
% There can be between one and three candidates, for every pitch
% from the 0th to 12th fret.
% --------------------------------------------------------
% [string, fret] = icassp19_obtain_pitch_candidates(observedPitch, frequencyReference)
% --------------------------------------------------------
function [stringCandidates, fretCandidates] = icassp19_obtain_pitch_candidates(observedPitch, frequencyReference)

% find the nearest note based on log frequency
[~, candidates(1,2)] = min(min(abs(log(frequencyReference)-log(observedPitch))));
[~, candidates(1,1)] = min(abs(log(frequencyReference(:,candidates(1,2)))-log(observedPitch)));

% translate to human understandable string and fret.
string = candidates(1,1);
candidates(1,2) = candidates(1,2)-1;
fret   = candidates(1,2);

% Find all candidates using the equal tempered scale
if string == 1 && (fret > 4 && fret < 10)
    candidates(2,1) = string+1;
    candidates(2,2) = fret-5;
end
if string == 1 && fret > 9
	candidates(2,1) = string+1;
    candidates(2,2) = fret-5;
    candidates(3,1) = string +2;
    candidates(3,2) = fret-10;
end
if string == 2 && fret < 5 
	candidates(2,1) = string-1;
	candidates(2,2) = fret+5;
end
if string == 2 && (fret > 4 && fret < 8)
	candidates(2,1) = string+1;
	candidates(2,2) = fret-5;
	candidates(3,1) = string-1;
	candidates(3,2) = fret+5;
end
if string == 2 && (fret > 7 && fret < 10)
	candidates(2,1) = string+1;
	candidates(2,2) = fret-5;
end
if string == 2 && (fret > 9)
	candidates(2,1) = string+1;
	candidates(2,2) = fret-5;
	candidates(3,1) = string+2;
	candidates(3,2) = fret-10;
end

if string == 3 && fret < 5 
	candidates(2,1) = string-1;
	candidates(2,2) = fret+5;
end
if string == 3 && fret < 3 
	candidates(3,1) = string-2;
	candidates(3,2) = fret+10;
end
if string == 3 && (fret > 4 && fret < 8)
	candidates(2,1) = string+1;
	candidates(2,2) = fret-5;
	candidates(3,1) = string-1;
	candidates(3,2) = fret+5;
end
if string == 3 && (fret > 7 && fret < 10)
	candidates(2,1) = string+1;
	candidates(2,2) = fret-5;
end
if string == 3 && (fret > 8)
	candidates(2,1) = string+1;
	candidates(2,2) = fret-5;
	candidates(3,1) = string+2;
	candidates(3,2) = fret-9;
end

if string == 4 && fret < 4 
	candidates(2,1) = string-1;
	candidates(2,2) = fret+5;
end
if string == 4 && fret < 3 
	candidates(3,1) = string-2;
	candidates(3,2) = fret+10;
end
if string == 4 && (fret > 3 && fret < 8)
	candidates(2,1) = string+1;
	candidates(2,2) = fret-4;
	candidates(3,1) = string-1;
	candidates(3,2) = fret+5;
end
if string == 4 && (fret > 7 && fret < 10)
	candidates(2,1) = string+1;
	candidates(2,2) = fret-4;
end
if string == 4 && (fret > 8)
	candidates(2,1) = string+1;
	candidates(2,2) = fret-4;
	candidates(3,1) = string+2;
	candidates(3,2) = fret-9;
end

if string == 5 && fret < 5 
	candidates(2,1) = string-1;
	candidates(2,2) = fret+4;
end
if string == 5 && fret < 4 
	candidates(3,1) = string-2;
	candidates(3,2) = fret+9;
end
if string == 5 && (fret > 4 && fret < 9)
	candidates(2,1) = string+1;
	candidates(2,2) = fret-5;
	candidates(3,1) = string-1;
	candidates(3,2) = fret+4;
end
if string == 5 && fret > 8
	candidates(2,1) = string+1;
	candidates(2,2) = fret-5;
end 

if string == 6 && fret < 4 
	candidates(2,1) = string-1;
	candidates(2,2) = fret+5;
	candidates(3,1) = string-2;
	candidates(3,2) = fret+9;
end
if string == 6 && (fret > 3 && fret < 8) 
	candidates(2,1) = string-1;
	candidates(2,2) = fret+5;
end
% translate to MATLAB indexing>0
candidates(:,2)=candidates(:,2)+1;
candidates=candidates';
% sort string-wise ascending
[~,ndx] = sort(candidates(1,:));
candidates = candidates(:,ndx);

	    fretCandidates   = candidates(2,:)-1;
	    stringCandidates = candidates(1,:);
	    numCandidates = length(stringCandidates);
end