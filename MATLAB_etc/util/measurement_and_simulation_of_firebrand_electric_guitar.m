clear all;
%close all;
for iiii=1:500
%% Initialise parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pluck position and displacement
P = 1 / 3; 	
force = .05;

L0 = 0.6411 + [1 2 9 4 5 2]'* 1e-3; % distance from nut to bridge

numStrings = 6;
strNdx = (1:numStrings)';
numFrets = 12;
fretNdx = (0:numFrets)';
    
%Material properties
Esteel = 2.27e11; % Young's modulus for steel
G = 79.3e9; % known shear modulus steel

rhoCore = 7950; %[kg/m^3]
rhoWrapping = 6000; %[kg/m^3]

dFull = [.010; .013; .0172; .026; .036; .046] * 0.0254; %full diameter for the  strings
dCore = [.010; .013; .0172; .0146; .016; .018] * 0.0254; %core diameter without wrapping

T0 = [16.5; 16; 17.5; 18.5; 20; 17;]*4.45;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add noise to length
L0 = L0.*2.^(-fretNdx'/12);
L0=L0+randn(size(L0)).*L0*0.005;

% Add noise to T0
T0=T0+randn(size(T0)).*T0*0.005;

% Add noise to dCore
dCore=dCore+randn(size(dCore)).*dCore*0.005;

% Add noise to dWrapping
dWrapping = (dFull - dCore) / 2; %diameter of wrapping wire
dWrapping=dWrapping+randn(size(dWrapping)).*dWrapping*0.005;

% Add noise to plucking force
force=force+randn(size(force)).*force*.005;




%%%%%%%%%%%%%%%%%%%%%
% calculate B and f0

ACore = (pi*(dCore/2).^2);% (pi * dCore.^2) / 4; %Crosssection core [m^2]

% Mass-per-unit length
mu = ACore .* rhoCore + rhoWrapping * ((2 * dWrapping + dCore).^2 - dCore.^2) * (pi / 4);

% Calculating deltaL from material properties
D = dCore + dWrapping;
TcOverTw = (8 * ACore .* D.^3 .* Esteel) ./ (G * dWrapping.^5);
Tc = T0 ./ ((1 ./ TcOverTw) + 1);
Tw = T0 ./ (TcOverTw + 1);
deltaL = (L0 .* Tc) ./ (ACore * Esteel + Tc);

%calculate effective Young's modulus for all strings
E = (Tc ./ ACore) ./ (deltaL ./ (L0 - deltaL));
Eeff = (T0 ./ ACore) ./ (deltaL ./ (L0 - deltaL));

% Transverse Displacement [Abbot, Strings in the 16th and 17th Centuries, 1974, Appendix 3]
deltaP = ((L0 * P .* force .* (1-P)) ./ T0).^2; 

% Length extension
deltaDeltaL = sqrt((P .* L0).^2 + deltaP.^2) + sqrt(((1 - P).*L0).^2 + deltaP.^2) - L0;

% Displacement as if plucked in the middle (retaining the length extension)
deltaHalf = sqrt((deltaDeltaL.^2 + deltaDeltaL .* L0 * 2) / 4);

% Calculate inharmonicity factor
K = (pi^3 * Eeff .* dCore.^2) ./ (16 * T0 .* L0.^2);

%K = -(pi^3 * Eeff .* dCore.^2) ./ (16 * T0 .* L0.^1);

% B is a sum of intrinsic and pluck deflection
BIntrinsic(:,:,iiii) = (K ./ 4) .* dCore.^2;
BPluck(:,:,iiii) = ((K * 3) ./ 8) .* deltaHalf.^2;

f0(:,:,iiii) = sqrt(T0./mu)./L0/2;

end
BMixed = BIntrinsic + BPluck;
BSimulations = BMixed;
pitchSimulations = f0;

w0 = f0;

%load ~/repositories/guitar_string_finger_and_pluck_estimation/util/inharmonicity/mats/BMean_Firebrand_40ms_wo.mat
load measurements.mat; %./mats/BMean_Firebrand_40ms_weightingSquared.mat
mf=mean(f0TableFirebrand,3); mb =mean(BTableFirebrand,3);

%
% 
% figure; hold on;   
% plot(mf(:),mb(:),'b*');        
% clear mu;
% mu1 = mean(pitchSimulations ,3);      ;  
% mu2 = mean(BSimulations,3);
% plot(mu1(:), mu2(:),'ko')  
% ylabel('B [m^{-2}]'); xlabel('pitch [Hz]');
% legend('$\widehat{\mu}_k$','Simulated ${\mu}_k$');
% title('Every \mu_k is computed from 10 observations');
% %%
% figure; plot((f0Table(:)),(BTable(:)),'k.'); grid minor; title('Estimates on Electric');
% hold on; plot(mf(:),mb(:),'b*')
% hold on; plot(median(f0Table,3),median(BTable,3),'r*');
% ylabel('B [m^{-2}]'); xlabel('pitch [Hz]');
% legend('\{$\widehat{\omega}_0$, $\widehat{B}$\}','$\mu_k$ (mean)','median');

%
% figure(20);  hold on;
% plot(w0(:),BMixed(:),'g.')
% xlabel('Simulation of pitch [Hz]')
% ylabel('Simulation of B ')
% title('Monte Carlo simulations')
% grid minor
%
figure(21); clf
scatter1 = scatter(w0(:),BMixed(:),'filled','.','linewidth',.1,'MarkerFaceColor',[0.43 0.43 0.43],'MarkerEdgeColor',[0.41 0.41 0.41]); 
% Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
scatter1.MarkerFaceAlpha = .5;
scatter1.MarkerEdgeAlpha = .5;
title('Feature space of the electric guitar')
hold on ; scatter((f0TableFirebrand(:)),(BTableFirebrand(:)),'.','MarkerFaceColor',[.9 .69 .0 ],'MarkerEdgeColor',[0.91 0.61 0.01]);

grid minor
ylabel('B [\cdot]'); xlabel('\omega_0 [Hz]');
legend('simulated','measured');
xlim([80 680]); ylim([1e-5 6.4e-4]);
%%
%figure(1); plot( squeeze(w0(1,:,:)),squeeze(BMixed(1,:,:)),'k.',squeeze(w0(2,:,:)),squeeze(BMixed(2,:,:)),'k.' ); grid minor; xlabel('Simulation of pitch [Hz]'); ylabel('Simulation of B '); title('Monte Carlo simulation')
%figure(2); plot( squeeze(w0(3,:,:)),squeeze(BMixed(3,:,:)),'k.',squeeze(w0(4,:,:)),squeeze(BMixed(4,:,:)),'k.' ); grid minor; xlabel('Simulation of pitch [Hz]'); ylabel('Simulation of B '); title('Monte Carlo simulation')
%figure(3); plot( squeeze(w0(5,:,:)),squeeze(BMixed(5,:,:)),'k.',squeeze(w0(6,:,:)),squeeze(BMixed(6,:,:)),'k.' ); grid minor; xlabel('Simulation of pitch [Hz]'); ylabel('Simulation of B '); title('Monte Carlo simulation')

