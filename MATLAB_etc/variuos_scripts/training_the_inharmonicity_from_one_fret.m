clear all; %close all;
fig_a = 0;
saveWorkspaceFile = 0;
addpath ~/repositories/guitar_string_finger_and_pluck_estimation/util/
mirverbose(0); mirwaitbar(0);
duration = 0.04; %sec
fs = 44100; N=round(fs*duration);
stringArea = 1:6;
fretArea = 6:6;
nFFT = 2^19;
guitarType = 'firebird';
betaRes = 1e-7;

% reomive after output
for fretArea=13:13

file = strcat('~/Dropbox/pluck_position/testfiles/',guitarType,'1','/string1','/','0.wav')
tmpInput = audioread(file);
tmpInput=tmpInput/max(abs(tmpInput));
% [value startidx] = smc_max(tmpInput.^2); 
% tmpInput = tmpInput(startidx+(round(fs*(0.002))):startidx+(round(fs*(0.002+duration))));
% sig = tmpInput(1:N)/max(abs(tmpInput(1:N)));

% read and segment at onset positions
[sig] = pluck_position_segment_from_all_onsets(tmpInput,fs,duration);

x = smc_hilbert(sig);
% define f0 for open strings
[f,Xf0] = smc_fft_LengthN(x,fs,nFFT);
openEStringF0 = smc_harmonic_summation_tuner_for_testing(Xf0, [35:350], 5, fs);
f0 = openEStringF0*[1 2^(5/12) 2^(10/12) 2^(15/12) 2^(19/12) 2^(24/12)];

numberOfFiles = 10;
%upliftFactor = 0.1;
%%  Read 6 strings 13 frets into "sig" matrix.
for numberOfFilesToRead = 1:numberOfFiles
    for string = stringArea
        for fret = fretArea
            file = strcat('~/Dropbox/pluck_position/testfiles/',guitarType,num2str(numberOfFilesToRead),'/string',num2str(string),'/',num2str(fret-1),'.wav');
            tmpInput = audioread(file);
            %[value startidx] = smc_max(tmpInput.^2); 
            %tmpInput = tmpInput(startidx+(round(fs*(0.002))):startidx+(round(fs*(0.002+duration))));
            %sig(:,string,fret) = tmpInput(1:N)/max(abs(tmpInput(1:N)));
            tmpInput=tmpInput/max(abs(tmpInput));
            [sig(:,string,fret)] = pluck_position_segment_from_all_onsets(tmpInput,fs,duration);
            %sig(:,string,fret) = smc_fadein(sig(:,string,fret),7,'log');
            %sig(:,string,fret) = smc_fadeout(sig(:,string,fret),5,'lin');
            x(:,string,fret) = smc_hilbert(sig(:,string,fret));
            % define f0 for open strings
            if fret == 1
               [f,Xf0] = smc_fft_LengthN(x(:,string,fret),fs,2^19);
               f0(string) = smc_harmonic_summation_tuner_for_testing(Xf0, [35:350], 5, fs);
            end
        end
    end

%%
    for string = stringArea 
        for fret = fretArea
            x(:,string,fret) = smc_apply_gaussian_window(x(:,string,fret));
            [f,X(:,string,fret)] = smc_fft_LengthN(x(:,string,fret),fs,nFFT);
            % heuristic noise removal
            %[X(:, string, fret), binMeanVector, diffVEctor, MedianVector, threshVector] = smc_equalize_harmonics(X(:,string,fret), fs, 2^((fret-1)/12)*f0(string));
        end
    end
    for string = stringArea
        for fret = fretArea
            % estimate pitch and inharmonicity coefficient.
            betaCoeffSearchArea = [1e-5:betaRes:6e-4];

            [ff0,Xf0] = smc_fft_LengthN(x(:,string,fret),fs,nFFT);
            f0HS = smc_harmonic_summation_tuner_for_testing(X(:,string,fret), [f0(string)*2^((fret-1)/12)-10:f0(string)*2^((fret-1)/12)+10], 3, fs);
            %L = 23+floor(2500/f0HS);
            L = min(54,floor(fs/2/f0HS));
            %L = 23+floor(2500/f0HS); % for f0HS = 700; L must max be 26 then, due to sampling rate.
            [est_f0(numberOfFilesToRead,string,fret) betaCoeff(numberOfFilesToRead,string,fret)] ...
            = icassp19_inharmonic_summation(X(:,string,fret), f0HS, L, fs, betaCoeffSearchArea,nFFT);

            %= smc_inharmonic_summation_tuner(X(:,string,fret),f0HS,L,fs,betaCoeffSearchArea, nFFT);
        end
    end
    % Visualize the betaCoeff for each string and fret
    for string=stringArea
        for fret=fretArea
            betaTable(1,string,numberOfFilesToRead)=(betaCoeff(1,string,fret));
            f0Table(1,string)=(est_f0(1,string,fret));
        end
    end
end
betaMean = mean(betaTable,3);

% build model
eta0 = fretArea-1;

%load /home/jmhh/repositories/guitar_string_finger_and_pluck_estimation/mats/betaMean_Firebird_40.mat;
%whos
fretIndex = [0 1 2 3 4 5 6 7 8 9 10 11 12]';
eta1 = fretIndex-eta0;
string = stringArea';%[1 2 3 4 5 6]';
 
betaModelApproximation = 2.^(eta1/6) * betaMean(1,string);
betaMean = betaModelApproximation;

figure(2), hold on;
    plot(fretIndex,betaModelApproximation');
    %xlim(fretIndex);
    xlabel('Fret'); ylabel('B'); 
    legend('str. 6','str. 5','str. 4','str. 3','str. 2','str. 1','model', 'Location', 'northwest'); 
    grid on;        
    %title('Western Gtr. with 0.012 inch.'); grid on;
    
if saveWorkspaceFile==1,
    outputFileName = strcat('~/repositories/guitar_string_finger_and_pluck_estimation/mats/betaMean',guitarType,'_40ms_from_',sprintf('%1.0fth_fret_betares%1.2fu_nFFT2^%1.0f',fretArea-1,betaRes*1e6,log2(nFFT)),'.mat')
    save(outputFileName,'betaTable','est_f0','betaMean');
    load(strcat('~/repositories/guitar_string_finger_and_pluck_estimation/mats/betaMean_',guitarType,'_40ms.mat'));
    figure(100)
        plot(betaMean)
        xlabel('Fret'); ylabel('B'); 
        legend('str. 6','str. 5','str. 4','str. 3','str. 2','str. 1','model', 'Location', 'northwest'); 
        grid on;  
    hold on;
        plot(betaModelApproximation,'k')

end
% for outputting B-files
sig=[];
x=[];
betaTable=[];
end