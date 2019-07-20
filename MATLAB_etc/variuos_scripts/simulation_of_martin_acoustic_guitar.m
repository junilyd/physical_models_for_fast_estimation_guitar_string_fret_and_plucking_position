clear all;
%close all;
tic;
for iiii=1:500
%% Initialise parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pluck position and displacement
P = 1 / 3; 	
force = .05;

L0 = 25.340*0.0254 ;%0.6411 + [1 2 9 4 5 2]'* 1e-3; % distance from nut to bridge

numStrings = 6;
strNdx = (1:numStrings)';
numFrets = 12;
fretNdx = (0:numFrets)';
    
%Material properties
Esteel = 2.27e11; % Young's modulus for steel
G = 79.3e9; % known shear modulus steel

rhoCore = 7950; %[kg/m^3]
rhoWrapping = 6250; %[kg/m^3]

%dFull = [.010; .013; .0172; .026; .036; .046] * 0.0254; %full diameter for the  strings
%dCore = [.010; .013; .0172; .0146; .016; .018] * 0.0254; %core diameter without wrapping

dFull = [.0115; .0151; .023 ; .032; .0445; .0564] * 0.0254; %full diameter for the  strings
dCore = [.0115; .0151; .0136; .014; .0153; .0185] * 0.0254; %core diameter without wrapping

T0 = ([12.8; 12; 17.5; 18.5; 20; 17;]+9)*4.45;

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
toc;
BMixed = BIntrinsic + BPluck;
BSimulations = BMixed;
pitchSimulations = f0;

w0 = f0;

load ~/repositories/guitar_string_finger_and_pluck_estimation/util/inharmonicity/mats/BMean_Martin_40ms_wo.mat% Firebrand_40ms_wo.mat
mf=mean(f0Table,3); mb =mean(BTable,3);

minX = 305; maxX = 336;
minY = .68e-4; maxY = 1.26e-4;
figure(22); clf

scatter1 = scatter(w0(:),BMixed(:),10,'filled','o','MarkerFaceColor',0*[0.43 0.43 0.43],'MarkerEdgeColor',0*[0.21 0.21 0.21]); 
% Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
%scatter1.MarkerFaceAlpha = .4;
%scatter1.MarkerEdgeAlpha = .5;
hold on ; scatter((f0Table(:)),(BTable(:)),6,'o','MarkerFaceColor',[.9 .69 .0 ],'MarkerEdgeColor',[0.91 0.61 0.01]);
%hold on; scatter(oneClassf0(:),oneClassB(:),'o','MarkerFaceColor',[.9 .69 .0 ],'MarkerEdgeColor',[0.91 0.61 0.01]);

plot([minX minX],[minY maxY],'k',[maxX maxX],[minY maxY],'k')
plot([minX maxX],[minY minY],'k',[minX maxX],[maxY maxY],'k')
plot([minX 200],[minY 2.98e-4],'--k')
plot([maxX 470],[minY 2.98e-4],'--k')
plot([minX 200],[maxY 4.3e-4],'--k')
plot([maxX 470],[maxY 4.3e-4],'--k')
l=legend('Simulated','Measured');
set(l,'Position',[.63 .8470 0.1 0.2],'orientation', 'horizontal'); legend boxoff %,,"orientation", "horizontal");legend boxoff 

%set(gca,'Color',0.95*[1 1 1])

grid minor
ylabel('B [\cdot]'); xlabel('\omega_0 [Hz]');
xlim([80 500]); ylim([1e-5 5e-4]);
% zoom in
axes('position',[.35 .6 .5 .25])

w0=w0(:);
BMixed = BMixed(:);
f0Table = f0Table(:);
BTable = BTable(:);
oneClassw0 = w0((w0>minX) & (w0<maxX) & BMixed>minY & BMixed<maxY);
oneClassB1 = BMixed(w0>minX & (w0<maxX) & BMixed>minY & BMixed<maxY);
oneClassf0 = f0Table(f0Table>minX & (f0Table<maxX) & BTable>minY & BTable<maxY);
oneClassB = BTable(f0Table>minX & (f0Table<maxX) & BTable>minY & BTable<maxY);

% scatter1 = scatter(oneClassw0(:),oneClassB1(:),'filled','o','MarkerFaceColor',[0.43 0.43 0.43],'MarkerEdgeColor',[0.41 0.41 0.41]) % plot on new axes
% scatter1.MarkerFaceAlpha = .3;
% scatter1.MarkerEdgeAlpha = .5;

%scatter(w0(:),BMixed(:),'filled','*','MarkerFaceColor',[0.43 0.43 0.43],'MarkerEdgeColor',[0.41 0.41 0.41]); 
scatter1 = scatter(oneClassw0(:),oneClassB1(:),'filled','o','MarkerFaceColor',0*[0.43 0.43 0.43],'MarkerEdgeColor',0*[0.21 0.21 0.21]) % plot on new axes
% scatter1.MarkerFaceAlpha = .3;
% scatter1.MarkerEdgeAlpha = .5;
hold on; scatter(oneClassf0(:),oneClassB(:),'o','MarkerFaceColor',[1 1 1 ],'MarkerEdgeColor',[1 1 1]);
axis tight; box on;% set(gca,'Color',0.92*[1 1 1])
hold on; scatter(oneClassf0(:),oneClassB(:),'o','MarkerFaceColor',1.1*[.9 .69 .0 ],'MarkerEdgeColor',[0.91 0.61 0.01]);

