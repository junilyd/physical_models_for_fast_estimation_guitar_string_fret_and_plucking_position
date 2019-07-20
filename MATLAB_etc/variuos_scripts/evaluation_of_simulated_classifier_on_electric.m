clear all;
addpath(genpath('../util'));
mirverbose(0);
addpath mats

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
L0=L0+0.005*randn(size(L0)).*L0;

% Add noise to T0
T0=T0+0.005*randn(size(T0)).*T0;

% Add noise to dCore
dCore=dCore+0.005*randn(size(dCore)).*dCore;

% Add noise to dWrapping
dWrapping = (dFull - dCore) / 2; %diameter of wrapping wire
dWrapping=dWrapping+0.005*randn(size(dWrapping)).*dWrapping;

% Add noise to plucking force
force=force+0.005*randn(size(force)).*force;


%%%%%%%%%%%%%%%%%%%%%
% calculate B and w0

 %Crosssection core [m^2]
ACore = (pi*(dCore/2).^2);

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

%% Calculate inharmonicity factor
K = (pi^3 * Eeff .* dCore.^2) ./ (16 * T0 .* L0.^2);

% B is a sum of intrinsic and pluck deflection
BIntrinsic(:,:,iiii) = (K ./ 4) .* dCore.^2;
BPluck(:,:,iiii) = ((K * 3) ./ 8) .* deltaHalf.^2;

f0(:,:,iiii) = sqrt(T0./mu)./L0/2;

end
BMixed = BIntrinsic + BPluck;
BSimulations = BMixed;
pitchSimulations = f0;
w0 = f0;

figure(28); hold on; plot((f0(:)),(BMixed(:)),'.');

% test the classifier
size(w0)
size(BMixed)
K = size(w0,1)*size(w0,2)
kk=0;
featureMatrix=[];
model=[];
for strings=6:-1:1
    for frets=1:13
    kk=kk+1;
    featureMatrix(:,:,kk) = [squeeze(w0(strings,frets,:)) squeeze(BMixed(strings,frets,:))];
    model.mu(:,kk) = mean(featureMatrix(:,:,kk));
    model.Sigma(:,:,kk) = cov(featureMatrix(:,:,kk));
    model.w(kk)=1/K;
    end
end

% capture an observation
load measurements.mat%   BMean_Firebrand_40ms.mat	% to be used as previously computed features (w0,B) using data and estimator of [25].
kk=0;
errCounts = 0;
for BB = 1:size(BTableFirebrand,3)
    w0 = f0TableFirebrand(:,:,BB)';
    B = BTableFirebrand(:,:,BB)';

    for trueString = 1:6,
        for trueFret = 0:12

        x = [w0(trueString,trueFret+1); B(trueString,trueFret+1)];

        [sC, fC] = icassp19_obtain_pitch_candidates(x(1),w0);
        ndx = (sC*13)-12+fC;
        cntCandidates=0;
        
        for ndx = (sC*13)-12+fC
            mu = model.mu(:,ndx);
            C = model.Sigma(:,:,ndx);
            P = model.w(ndx);

            %% Classifier
            cntCandidates=cntCandidates+1;
            [J(cntCandidates)] = max(-log(det(C)) + 2*log(P) ...
                               - mu'*inv(C)*mu + 2*x'*inv(C)*mu - x'*inv(C)*x);
%                           
        end
        [~,I] = max(J); 
        clear J euclideanDistance;
        
        kk=kk+1;
        stringEstimate(kk) = sC(I);
        fretEstimate(kk) = fC(I);

        if fretEstimate(kk) == -1, fretOptions=[1:13]; fretEstimate(kk)=fretOptions(end-1); stringEstimate(kk)=stringEstimate(kk)-1;end
            if trueFret~=fretEstimate(kk)
                errCounts = errCounts+1;
                fprintf('--------------\n%5.0f %5.0f\n%5.0f %5.0f\n--------------\n' ...
                    ,trueString,stringEstimate(kk),trueFret,fretEstimate(kk))
            end
        end
    end
    figure(28); plot((w0),(B),'k.'), hold on;
end
    xlabel('\omega_0  [Hz]'); 
    ylabel('$B [\cdot]$','Interpreter','latex');