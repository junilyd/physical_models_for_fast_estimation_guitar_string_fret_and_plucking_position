function [stringEstimate, fretEstimate] = evaluation_of_simulated_classifier_on_electric_function(x)
%clear all;
addpath ~/repositories/guitar_string_finger_and_pluck_estimation/icassp_plucking_estimation/MATLAB/util/
%close all;
for iiii=1:5000
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
BMixed = BIntrinsic + BPluck;
BSimulations = BMixed;
pitchSimulations = f0;

w0 = f0;


%figure(28); hold on; plot((f0(:)),(BMixed(:)),'.');

% test the classifier
size(w0);
size(BMixed);
K = size(w0,1)*size(w0,2);
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
%load ~/repositories/guitar_string_finger_and_pluck_estimation/util/inharmonicity/mats/BMean_Firebrand_40ms_weightingSquared.mat	

%kk=0;
%errCounts = 0;
%for BB = 1:size(BTable,3)
%    w0 = f0Table(:,:,BB)';
%    B = BTable(:,:,BB)';

    %for trueString = 1:6,
    %    for trueFret = 0:12

        %x = [w0(trueString,trueFret+1); B(trueString,trueFret+1)];
        % x should be a vector with w0 and B
       
        w0Reference = vec2mat(model.mu(1,:),13);
        [sC, fC] = icassp19_obtain_pitch_candidates(x(1),w0Reference);
        %ndx = (sC*13)-12+fC;
        cntCandidates=0;
        
        % x = x./[675 6.5e-4]'; % for Euclidean distance 

        for ndx = (sC*13)-12+fC
            mu = model.mu(:,ndx);
            %mu = mu./[675 6.5e-4]'; % for euclidean distance
            C = model.Sigma(:,:,ndx);
            P = model.w(ndx);

            %% Classifier
            cntCandidates=cntCandidates+1;
            [J(cntCandidates)] = max(-log(det(C)) + 2*log(P) ...
                               - mu'*inv(C)*mu + 2*x'*inv(C)*mu - x'*inv(C)*x);
%                           
%             euclideanDistance(cntCandidates)   = norm(x-mu);

        end
        [~,I] = max(J); 
%        [~,I] = min(euclideanDistance); % for Euclidean distance

        clear J euclideanDistance;
%             stringEstimate(trueString,fretNdx,recDir) = sC(I);
%             fretEstimate(trueString,fretNdx,recDir) = fC(I);
        
        
        %kk=kk+1;
        stringEstimate = sC(I);
        fretEstimate = fC(I);

%         if fretEstimate(kk) == -1, fretOptions=[1:13]; fretEstimate(kk)=fretOptions(end-1); stringEstimate(kk)=stringEstimate(kk)-1;end
%             if trueFret~=fretEstimate(kk)
%                 errCounts = errCounts+1;
%                 fprintf('--------------\n%5.0f %5.0f\n%5.0f %5.0f\n--------------\n' ...
%                     ,trueString,stringEstimate(kk),trueFret,fretEstimate(kk))
%             end
%         end
    %end
%    figure(28); plot((w0),(B),'k.'), hold on;
%end
%    xlabel('\omega_0  [Hz]'); 
%    ylabel('$B [\cdot]$','Interpreter','latex');