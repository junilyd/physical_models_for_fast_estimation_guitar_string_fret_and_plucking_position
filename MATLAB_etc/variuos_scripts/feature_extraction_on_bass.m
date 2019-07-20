clear all; %close all;
fig_a = 0; % plot ? yes: 1
saveWorkspaceFile = 1; % save *.mat file with model(s) ? yes: 1

stringOptions = 1:4;
fretOptions = 1:19;

addpath ~/repositories/guitar_string_finger_and_pluck_estimation/util/
addpath(genpath('util'));
mirverbose(0); mirwaitbar(0);
SNROptions = [180 90 80 70 60 50 40 30 20 10 0]%[80 50 30 10 5 0 -5 -10 -15 -20];
guitarType = 'fenderBass'%'firebird';% 'martin';

% Initialize implementation constants 
segmentDuration = 40e-3; % segment duration in seconds.
fs = 44100; % N=round(fs*segmentDuration);
M = 25; % assumed high number of harmonics (M>>1).
MInitial = 3; % number of harmonics for initial harmonic estimate (B=0).
f0Limits = [75 700]/2; % boundaries for f0 search grid in Hz.
nFFT = 2^19; % Length of  zero-padded FFT.
betaRes = 1e-7; % resolution of search grid for B in m^(-2)
BSearchGrid = [1e-5:betaRes:6.6e-4]; % searh grid for B. 

%%
%for SNR = SNROptions
for fretCounter=fretOptions
%% Read 6 strings 13 frets into "sig" matrix.
    % Feature Extraction (pitch and inharmonicity)
for numberOfFilesToRead = 1:10
    for string = stringOptions
        for fret = fretCounter
            % read audio
            recordingPath = strcat('~/Dropbox/pluck_position/testfiles/',guitarType,num2str(numberOfFilesToRead),'/string',num2str(string),'/',num2str(fret-1),'.wav');
            [recording,fs] = audioread(recordingPath);
            recording=recording/max(abs(recording));
            [sig] = icassp19_segment_from_all_onsets(recording,fs,segmentDuration);
            sig = sig(:,1);
            %sig = smc_AWGN(sig, SNR);
            % Hilbert transform and windowing
            x(:,string,fret) = icassp19_hilbert_transform(sig,fs);
            x(:,string,fret) = icassp19_apply_gaussian_window(x(:,string,fret));
            %% Feature extraction with the inharmonic pitch estimator
            % The implementation of Eq. (17) is done with one FFT, since it is fast
            % Hence, it is equivalent to harmonic summation which the in the proposed 
            % method is extended to inharmonic summation.
            % See details on harmonic summation in Christensen and Jakobsson [27].
            [~,X(:,string,fret)] = icassp19_fft(x(:,string,fret), fs, nFFT);
            omega0Initial = icassp19_harmonic_summation(X(:,string,fret), f0Limits, MInitial, fs);
            %M = min(54,floor(fs/2/omega0Initial));
            [pitchEstimate(numberOfFilesToRead,string,fret), BEstimate(numberOfFilesToRead,string,fret)] ...
            = icassp19_inharmonic_summation(X(:,string,fret), omega0Initial, M, fs, BSearchGrid,nFFT)
        end
    end
end

%% training the model.
for string = stringOptions 
    for fret = fretCounter
        % compute expected values
        expectedPitch(1,string,fret) = mean(pitchEstimate(:,string,fret));
        expectedB(1,string,fret) = mean(BEstimate(:,string,fret));
    end
end
for string = stringOptions 
    for fret = fretCounter
        % build a model in original units (Hz and m^(-2))
        % generalized fret settings
        f0 = fretCounter-1;
        fretIndex = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18]';
        f1 = fretIndex-f0;
        stringNdx = stringOptions';
        % model in original units
        BModelOriginalUnits(:,:,fret) = 2.^(f1/6) * expectedB(1,stringNdx,fret);
        pitchModelOriginalUnits(:,:,fret) = 2.^(f1/12) * expectedPitch(1,stringNdx,fret);
        % compute normalization factor
        normalizationConst(1,fret) = max(max(abs(pitchModelOriginalUnits(:,:,fret))));
        normalizationConst(2,fret) = max(max(abs(BModelOriginalUnits(:,:,fret))));
        % compute mu and lambda
        muPitch(:,:,fret) = pitchModelOriginalUnits(:,:,fret)./normalizationConst(1,fret)
        muB(:,:,fret) = BModelOriginalUnits(:,:,fret)./normalizationConst(2,fret)  
        lambda(:,:,string,fret) = cov([muPitch(:,string,fret) muB(:,string,fret)])

    end   
end
%% define output variables for the model of a given fret.
for string = stringOptions 
    for fret = fretCounter      
        tmp1 = muPitch(:,:,fret);
        tmp2 = muB(:,:,fret) ;
        tmp3(string) = mean(diag(lambda(:,:,string,fret)));
        mu = [tmp1(:) tmp2(:)];
        sigma = mean(tmp3);
        % labels ?
        labels = [(1:length(mu))'];
        fretLabel = repmat([1:fretOptions(end)]'-1,length(stringOptions),1);
        stringLabel = reshape([stringOptions].*ones(fretOptions(end),length(stringOptions)),length(stringOptions)*fretOptions(end),1);
        labels = [labels, stringLabel, fretLabel];
        normalizationConstant = [normalizationConst(1,fret) normalizationConst(2,fret)];
    end
end

%%
if saveWorkspaceFile==1,
for string = stringOptions 
    for fret = fretCounter
        clear tmp1 tmp2 x sig f1 fig_a fret string fretIndex omega0Initial recording stringNdx X labels
        outputFileName = strcat('/home/jmhh/WASPAA2019/mats/trained_model_of_',sprintf(guitarType),'_from_',sprintf('%1.0fth_fret',fretCounter-1),'.mat')
        save(outputFileName);
    end
    end
end

%end
end