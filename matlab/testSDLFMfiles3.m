%clc
clear
% randn('state', 1e6)
% rand('twister', 1e6)
randn('state', 1967688546)
rand('twister', 1967688546)
%numSeed = rand('twister');
%numSeed = numSeed(1);
% Create sdlfmKern1

nIntervals = 4;
tInit = 0;
tFinal = 5;
nPoints = 50;
nPoints2 = 40;
t1 = linspace(tInit,tFinal, nPoints)';
t2 = linspace(tInit,tFinal, nPoints2)';
options.nIntervals = nIntervals;
options.switchingTimes = [-0.1 2 1 0.5];
%options.switchingTimes = [-0.1  2];
options.nlfPerInt = 2;

LIC = [1 0;0 1];
covIC = [1 0;0 1];

sdlfmKern1 = kernCreate(t1, {'parametric', options ,'sdlfm'});
[params, names] = kernExtractParam(sdlfmKern1);
inverseWidthIndx = paramNameRegularExpressionLookup(sdlfmKern1, 'inverse width', 1);
inverseWidth = 1 + 2*rand(1, length(inverseWidthIndx));
sensitivityIndx = paramNameRegularExpressionLookup(sdlfmKern1, 'sensitivity', 1);
sensitivity = 3+ 5*rand(1,length(sensitivityIndx));
params(sensitivityIndx) = sensitivity; 
sdlfmKern1 = kernExpandParam(sdlfmKern1, params);
sdlfmKern1.inverseWidth = reshape(inverseWidth, options.nlfPerInt, options.nIntervals);
sdlfmKern1.mass = 0.1;
sdlfmKern1.spring = 0.5;
sdlfmKern1.damper = 1;

params = sdlfmKernExtractParam(sdlfmKern1);
sdlfmKern1 = sdlfmKernExpandParam(sdlfmKern1, params);

% Create sdlfmKern2
sdlfmKern2 = kernCreate(t1, {'parametric', options ,'sdlfm'});
sdlfmKern2.inverseWidth = sdlfmKern1.inverseWidth;
sdlfmKern2.sensitivity = 2*rand(size(sdlfmKern2.sensitivity));
sdlfmKern2.mass = 0.1;
sdlfmKern2.spring = 0.2;
sdlfmKern2.damper = 1;

params = sdlfmKernExtractParam(sdlfmKern2);
sdlfmKern2 = sdlfmKernExpandParam(sdlfmKern2, params);

t1 = t2;
sdlfmKern1 = sdlfmKern2;


covGrad = ones(length(t1), length(t2));


fhandle1 = str2func('sdlfmXsdlfmKernCompute');
fhandle2 = str2func('sdlfmXsdlfmKernGradient');

K = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% A detour to test gradients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PARAMETERS SYSTEM 1
%%%%% MASS 1
epsilon = 1e-6;
mass = sdlfmKern1.mass;
sdlfmKern1.mass = mass + epsilon;
K1 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);
sdlfmKern1.mass = mass - epsilon;
K2 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);
mass1n = 0.5*sum(sum((K1 - K2)))/epsilon;
sdlfmKern1.mass = mass;
%%%%%%%%%%% Helps to see where is a mistake
dim1 = zeros(1, options.nIntervals);
spVector = [cumsum(sdlfmKern1.switchingTimes) t1(end)+50];
for i =1:options.nIntervals
    newt1 = t1(t1> spVector(i) & t1<spVector(i+1));
    dim1(i) = length(newt1);
end
gradNum = zeros(length(dim1));
start1 = 1;
end1 =0;
for i=1:length(dim1)
    end1 = end1 + dim1(i);
    start2 = 1;
    end2 = 0;
    for j=1:length(dim1)
        end2 = end2 + dim1(j);
        gradNum(i,j) = 0.5*(sum(sum(K1(start1:end1, start2:end2) - K2(start1:end1, start2:end2))))/epsilon;
        start2 = end2 + 1;
    end
    start1 = end1 + 1;
end
%%%% SPRING 1
spring = sdlfmKern1.spring;
sdlfmKern1.spring = spring + epsilon;
K1 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);
sdlfmKern1.spring = spring - epsilon;
K2 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);
spring1n = 0.5*sum(sum((K1 - K2)))/epsilon;
sdlfmKern1.spring = spring;
%%%% DAMPER 1
damper = sdlfmKern1.damper;
sdlfmKern1.damper = damper + epsilon;
K1 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);
sdlfmKern1.damper = damper - epsilon;
K2 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);
damper1n = 0.5*sum(sum((K1 - K2)))/epsilon;
sdlfmKern1.damper = damper;
%%%% INVERSE WIDTH 1
indexIW = 4;
invW = sdlfmKern1.inverseWidth(indexIW);
sdlfmKern1.inverseWidth(indexIW) = invW + epsilon;
sdlfmKern2.inverseWidth(indexIW) = invW + epsilon;
K1 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);
sdlfmKern1.inverseWidth(indexIW) = invW - epsilon;
sdlfmKern2.inverseWidth(indexIW) = invW - epsilon;
K2 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);
inv1n = 0.5*sum(sum((K1 - K2)))/epsilon;
sdlfmKern1.inverseWidth(indexIW) = invW;
sdlfmKern2.inverseWidth(indexIW) = invW;
%%%% SENSITIVITY 1
indexS = indexIW;
senS = sdlfmKern1.sensitivity(indexS);
sdlfmKern1.sensitivity(indexS) =  senS + epsilon;
K1 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);
sdlfmKern1.sensitivity(indexS) =  senS - epsilon;
K2 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);
sensn = 0.5*(sum(sum(K1 - K2)))/epsilon;
sdlfmKern1.sensitivity(indexS) =  senS;
%%%% SWITCHING TIME 
indexSP = indexIW;
spTimes = sdlfmKern1.switchingTimes(indexSP);
sdlfmKern1.switchingTimes(indexSP) = spTimes + epsilon;
sdlfmKern2.switchingTimes(indexSP) = spTimes + epsilon;
K1 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);
sdlfmKern1.switchingTimes(indexSP) = spTimes - epsilon;
sdlfmKern2.switchingTimes(indexSP) = spTimes - epsilon;
K2 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);
sp1n = 0.5*sum(sum((K1 - K2)))/epsilon;
sdlfmKern1.switchingTimes(indexSP) = spTimes;
sdlfmKern2.switchingTimes(indexSP) = spTimes;
%%%%% PARAMETERS SYSTEM 2
%%%%% MASS 2
epsilon = 1e-6;
mass = sdlfmKern2.mass;
sdlfmKern2.mass = mass + epsilon;
K1 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);
sdlfmKern2.mass = mass - epsilon;
K2 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);
mass2n = 0.5*sum(sum((K1 - K2)))/epsilon;
sdlfmKern2.mass = mass;
%%%%%%%%%%% Helps to see where is a mistake
% dim1 = zeros(1, options.nIntervals);
% for i =1:options.nIntervals
%     newt1 = t1(t1> spVector(i) & t1<spVector(i+1));
%     dim1(i) = length(newt1);
% end
% gradNum = zeros(length(dim1));
% start1 = 1;
% end1 =0;
% for i=1:length(dim1)
%     end1 = end1 + dim1(i);
%     start2 = 1;
%     end2 = 0;
%     for j=1:length(dim1)
%         end2 = end2 + dim1(j);
%         gradNum(i,j) = 0.5*(sum(sum(K1(start1:end1, start2:end2) - K2(start1:end1, start2:end2))))/epsilon;
%         start2 = end2 + 1;
%     end
%     start1 = end1 + 1;
% end
%%%% SPRING 2
spring = sdlfmKern2.spring;
sdlfmKern2.spring = spring + epsilon;
K1 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);
sdlfmKern2.spring = spring - epsilon;
K2 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);
spring2n = 0.5*sum(sum((K1 - K2)))/epsilon;
sdlfmKern2.spring = spring;
%%%% DAMPER 2
damper = sdlfmKern2.damper;
sdlfmKern2.damper = damper + epsilon;
K1 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);
sdlfmKern2.damper = damper - epsilon;
K2 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);
damper2n = 0.5*sum(sum((K1 - K2)))/epsilon;
sdlfmKern2.damper = damper;
%%%% SENSITIVITY 2
indexS = indexIW;
senS = sdlfmKern2.sensitivity(indexS);
sdlfmKern2.sensitivity(indexS) =  senS + epsilon;
K1 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);
sdlfmKern2.sensitivity(indexS) =  senS - epsilon;
K2 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);
sens2n = 0.5*(sum(sum(K1 - K2)))/epsilon;
sdlfmKern2.sensitivity(indexS) =  senS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g1, g2, covGradLocal] = fhandle2(sdlfmKern1, sdlfmKern2, t1, t2, covGrad, covIC);

indexi = 2;
indexj = 2;
entry = LIC(indexi,indexj);
LIC(indexi,indexj) = entry + epsilon;
covIC = LIC*LIC';
%entry = covIC(indexi,indexj);
%covIC(indexi,indexj) = entry + epsilon;
K1 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);
LIC(indexi,indexj) = entry - epsilon;
covIC = LIC*LIC';
%covIC(indexi,indexj) = entry - epsilon;
K2 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);
licn = 0.5*(sum(sum(K1-K2)))/epsilon;
LIC(indexi,indexj) = entry;
covIC = LIC*LIC';
%covIC(indexi,indexj) = entry;
%%%%%%%%%%%%
% dim1= [20 15 15];
% gradNum = zeros(length(dim1));
% start1 = 1;
% end1 =0;
% for i=1:length(dim1)
%     end1 = end1 + dim1(i);
%     start2 = 1;
%     end2 = 0;
%     for j=1:length(dim1)
%         end2 = end2 + dim1(j);
%         gradNum(i,j) = 0.5*(sum(sum(K1(start1:end1, start2:end2) - K2(start1:end1, start2:end2))))/epsilon;
%         start2 = end2 + 1;
%     end
%     start1 = end1 + 1;
% end
%%%%%%%%%%%%%

[g1, g2, covGradLocal] = fhandle2(sdlfmKern1, sdlfmKern2, t1, t2, covGrad, covIC);


% Jij = zeros(2);
% Jji = zeros(2);
% gIC = zeros(2);
% for i=1:2
%     for j=1:2
%         Jij(i,j) = 1; Jji(j,i) = 1;
%         gIC(i,j) = sum(sum((LIC*Jij + Jji*LIC').*covGradLocal));
%         Jij(i,j) = 0; Jji(j,i) = 0;       
%     end
% end


%K2 = fhandle1(sdlfmKern1, sdlfmKern2, t1, t2, covIC);

% y = gaussSamp(K2, 5);
% 
% close all
% plot(t1, y);
% propa = get(gca);
% 
% for j=2:nIntervals
%     hold on
%     plot([sdlfmKern1.options.stimes(j) sdlfmKern1.options.stimes(j)], ...
%         propa.YLim, 'k--');    
% end
% 
% fileName = ['sdlfmSample' num2str(numSeed)];
% print('-dpdf', ['./figures/' fileName]);
% print('-depsc', ['./figures/' fileName]);





