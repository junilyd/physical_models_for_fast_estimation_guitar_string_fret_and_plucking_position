%% Function to count errors in estimated vs true
% v 1.0 (150324)
%
% INPUT 
%       errorCounter : total number errors
%        trueString  : true string in signal
%          trueFret  : true fret in signal
%  estimatedString   : estimated string
%    estimatedFret   : estimated fret
%       errorMatrix  : Matrix with location of errors
%
% OUTPUT
%        errorCounter: total number of errors
%         errorMatrix: Matrix with location of errors
%       
function [errorMatrix, errorCounter] = icassp19_count_and_localize_errors(trueString, trueFret, estimatedString, estimatedFret, errorMatrix, errorCounter)

        if trueString ~= estimatedString || trueFret ~=(estimatedFret)
			errorCounter = errorCounter+1;
			errorMatrix(trueString,trueFret+1) = errorMatrix(trueString,trueFret+1) + 1;
        else
            errorCounter=errorCounter;
        end
        

end