%fileNameMartin = 'mats/results_BFretModel_w_plausibilty_SNR0dB_no_weightingmartin_40ms_from_12th_fret_betares0.10u_nFFT2^19.mat';
fileNameMartin = 'mats/results_BFretModel_w_plausibilty_SNR0dB_no_weightingmartin_40ms_from_12th_fret_betares0.10u_nFFT2^19.mat';
fileNameFirebrand = 'mats/results_BFretModel_w_plausibilty_SNR0dB_no_weightingfirebird_40ms_from_12th_fret_betares0.10u_nFFT2^19.mat';

load(fileNameFirebrand)


errorRateTrainedElectric = errorRate;
errorRateSimElectric = errorRateSim;

load(fileNameMartin)
figure(100); 
plot(SNROptions, errorRateTrainedElectric, '-.ko', SNROptions, errorRateSimElectric, '-.k.', ...
    SNROptions, errorRate+eps, '-ko', SNROptions, errorRateSim+eps, '-k.');
xlabel('SNR [dB]');; ylabel('Error Rate');
legend('Trained model (electric)', 'Simulated Model (electric)','Trained model (acoustic)', 'Simulated Model (acoustic)','location','northwest')
set(gca,'xdir','reverse');
grid on;
xlim([0 50]);ylim([0 .12]);
