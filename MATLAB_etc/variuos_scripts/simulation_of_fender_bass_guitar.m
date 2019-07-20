clear all;
%load '/home/jmhh/WASPAA2019/mats/trained_model_of__from_12th_fret.mat';
%figure; plot(pitchEstimate(:,:,13),BEstimate(:,:,13)); hold on, plot(pitchModelOriginalUnits(:),BModelOriginalUnits(:),'.k')
%close all;
for iiii=1:500
%% Initialise parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pluck position and displacement
P = 1 / 3; 	
force = .05;

L0 = 34*0.0254 ;%0.6411 + [1 2 9 4 5 2]'* 1e-3; % distance from nut to bridge

numStrings = 4;
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

dFull = [.045;  .065 ;.085; .105;   ] * 0.0254; %full diameter for the  strings
dCore = [.019; .02; .022; .026; ] * 0.0254; %core diameter without wrapping
% dFull = [.0115; .0151; .023 ; .032; ] * 0.0254; %full diameter for the  strings
% dCore = [.0115; .0151; .0136; .014; ] * 0.0254; %core diameter without wrapping

T0 = [40.6; 49.0; 48.1; 42.4;]*4.45; %([12.8; 12; 17.5; 18.5;]+9)*4.45;

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

%load ~/repositories/guitar_string_finger_and_pluck_estimation/util/inharmonicity/mats/BMean_Martin_40ms_wo.mat% Firebrand_40ms_wo.mat
BModelOriginalUnits = [...
    0.1543    0.0674    0.0433    0.0451; ...
    0.1732    0.0756    0.0486    0.0506; ...
    0.1945    0.0849    0.0545    0.0568; ...
    0.2183    0.0953    0.0612    0.0638; ...
    0.2450    0.1069    0.0687    0.0716; ...
    0.2750    0.1200    0.0771    0.0804; ...
    0.3087    0.1347    0.0865    0.0902; ...
    0.3465    0.1512    0.0971    0.1013; ...
    0.3889    0.1697    0.1090    0.1137; ...
    0.4365    0.1905    0.1223    0.1276; ...
    0.4900    0.2138    0.1373    0.1432; ...
    0.5500    0.2400    0.1541    0.1608; ...
    0.6174    0.2694    0.1730    0.1805; ]*1.0e-03;
pitchModelOriginalUnits = [...
   41.2403   55.5997   74.1089   98.8028; ...
   43.6926   58.9059   78.5156  104.6779; ...
   46.2907   62.4086   83.1844  110.9023; ...
   49.0432   66.1196   88.1308  117.4969; ...
   51.9595   70.0513   93.3714  124.4837; ...
   55.0492   74.2168   98.9235  131.8859; ...
   58.3226   78.6299  104.8058  139.7282; ...
   61.7906   83.3055  111.0379  148.0369; ...
   65.4649   88.2591  117.6405  156.8396; ...
   69.3576   93.5073  124.6358  166.1658; ...
   73.4818   99.0675  132.0470  176.0465; ...
   77.8513  104.9583  139.8990  186.5148; ...
   82.4806  111.1995  148.2178  197.6055; ];

mf= pitchModelOriginalUnits(:);%  mean(f0Table,3); mb =mean(BTable,3);
mb = BModelOriginalUnits(:);

minX = 100; maxX = 115;
minY = .75e-4; maxY = 1.16e-4;
figure(22); clf

scatter1 = scatter(w0(:),BMixed(:),10,'filled','o','MarkerFaceColor',[0.43 0.43 0.43],'MarkerEdgeColor',[0.41 0.41 0.41]); 
% Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
scatter1.MarkerFaceAlpha = .4;
scatter1.MarkerEdgeAlpha = .5;
%%
hold on ; scatter((mf(:)),(mb(:)),20,'o','MarkerFaceColor',[.9 .69 .0 ],'MarkerEdgeColor',[0.91 0.61 0.01]);
plot([minX minX],[minY maxY],'k',[maxX maxX],[minY maxY],'k')
plot([minX maxX],[minY minY],'k',[minX maxX],[maxY maxY],'k')
plot([minX 98],[minY 2.98e-4],'--k')
plot([maxX 230],[minY 2.98e-4],'--k')
plot([minX 98],[maxY 4.3e-4],'--k')
plot([maxX 230],[maxY 4.3e-4],'--k')
l=legend('Simulated','Measured');
set(l,'Position',[.63 .8470 0.1 0.2],'orientation', 'horizontal'); legend boxoff %,,"orientation", "horizontal");legend boxoff 

grid minor
ylabel('B [\cdot]'); xlabel('\omega_0 [Hz]');
xlim([35 250]); ylim([1e-5 6e-4]);
hold on; plot(pitchModelOriginalUnits(:),BModelOriginalUnits(:),'.')

%% zoom in
axes('position',[.35 .6 .5 .25])
%
w0=w0(:);
BMixed = BMixed(:);
f0Table = mf(:);
BTable = mb(:);
oneClassw0 = w0((w0>minX) & (w0<maxX) & BMixed>minY & BMixed<maxY);
oneClassB1 = BMixed(w0>minX & (w0<maxX) & BMixed>minY & BMixed<maxY);
oneClassf0 = f0Table(f0Table>minX & (f0Table<maxX) & BTable>minY & BTable<maxY);
oneClassB = BTable(f0Table>minX & (f0Table<maxX) & BTable>minY & BTable<maxY);

%scatter(w0(:),BMixed(:),'filled','*','MarkerFaceColor',[0.43 0.43 0.43],'MarkerEdgeColor',[0.41 0.41 0.41]); 
scatter1 = scatter(oneClassw0(:),oneClassB1(:),'filled','o','MarkerFaceColor',[0.43 0.43 0.43],'MarkerEdgeColor',[0.41 0.41 0.41]) % plot on new axes
scatter1.MarkerFaceAlpha = .3;
scatter1.MarkerEdgeAlpha = .5;
hold on; scatter(oneClassf0(:),oneClassB(:),'o','MarkerFaceColor',[.9 .69 .0 ],'MarkerEdgeColor',[0.91 0.61 0.01]);
axis tight; box on;