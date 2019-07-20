% ----------------------------------------------------------------------
% Estimation of inharmonicity and pithc with inharmonic summation (approx. NLS)
%
%   INPUTS:
%           X:            fft of input signal recording
%           pitchInital : Initial estimate of pitch (usually with the harmonic assumption)
%           M:            Number of Harmonics to use.
%           fs:           Sample rate to use
%           BSearchGrid : Search grid for inharmonicity coefficient 
%   OUTPUTS:
%           pitch       : The pitch estimate
%           BEstimate   : The inharmonicity coefficient estimate
% ------------------------------------------------------------------------------------------------
% [pitchEstimate, BEstimate] = icassp19_inharmonic_summation(X, f0_area, L, fs,betaArea)
% ------------------------------------------------------------------------------------------------
function [pitchEstimate, BEstimate, costFunctionMaxVal] = icassp19_inharmonic_summation(X, pitchInitial, M, fs, BSearchGrid, nFFT)

pitchInitial = [pitchInitial 2*pitchInitial];
pitchWidth = (3+3)*nFFT/2^18*fs/nFFT; % (3+3) is for to cover the gaussian window mainlobe width

for pp = 1:1 % testing lower(1) or also upper(2) octave
    B=1;
    for BGrid=BSearchGrid
        pitchGrid = [pitchInitial(pp)-pitchWidth:fs/nFFT:pitchInitial(pp)+pitchWidth];
        for pII=1:length(pitchGrid)
            [HL] = icassp19_inharmonic_index(X, fs, M, pitchGrid(pII),BGrid);
            HL(HL>=length(X))=[];
            costFunction(:,B,pII) = sum(X(HL));
        end
        B=B+1;
    end
    
    betaGridSize  = size(costFunction,2);
    pitchGridSize = size(costFunction,3);
    %% CostFunction
    [C(pp) I] = max(costFunction(:));
    
    %% Pitch
    pitchIndex = ceil(I/betaGridSize);
    inharmonicPitch(pp) = pitchGrid(pitchIndex);

    %% B
    BGridIndex = mod(I,betaGridSize);
    if BGridIndex == 0, 
        BGridIndex = betaGridSize; 
    end
    BCandidates(pp) = BSearchGrid(BGridIndex);
end
%% loop end
[costFunctionMaxVal argMax] = max(C);
pitchEstimate = inharmonicPitch(argMax);
BEstimate = BCandidates(argMax);
end
