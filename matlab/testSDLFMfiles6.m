% TEST ALL THE DERIVATIVES FOR POSPOS, VELPOS, POSVEL, VELVEL, ACCELVEL,
% VELACCEL, ACCELACCEL

clc
clear
%randn('state', 1e6)
%rand('twister', 1e6)
numSeed = rand('twister');
numSeed = numSeed(1);
% Create sdlfmKern1
nIntervals = 3;
tInit = 0;
tFinal = 5;
nPoints = 50;

t1 = linspace(tInit,tFinal, nPoints)';


options.nIntervals = nIntervals;
options.switchingTimes = [-0.1  1  3];
options.nlfPerInt = 2;
options.isNormalised = 1;

%inverseWidth = 2*rand(1, options.nlfPerInt*options.nIntervals);

inverseWidth = ones(1, options.nlfPerInt*options.nIntervals);


% Create sdrbfKern

sdrbfKern = kernCreate(t1, {'parametric', options, 'sdrbf'});
[params, names] = sdrbfKernExtractParam(sdrbfKern);
indx = paramNameRegularExpressionLookup(sdrbfKern, 'inverse width .*', 1);
params(indx) = inverseWidth;
sdrbfKern = sdrbfKernExpandParam(sdrbfKern, params);

covGrad = ones(length(t1), length(t1));

K = sdrbfKernCompute(sdrbfKern, t1, t1);

%%%%% Test the gradients 
%%% INVERSE WIDTHS
epsilon = 1e-6;
iq = 2;
indexIW = 4;
invW = sdrbfKern.inverseWidth(indexIW);
sdrbfKern.inverseWidth(indexIW) = invW + epsilon;
K1 = sdrbfKernCompute(sdrbfKern, t1, t1);
sdrbfKern.inverseWidth(indexIW) = invW - epsilon;
K2 = sdrbfKernCompute(sdrbfKern, t1, t1);
inv1n = 0.5*sum(sum(K1{iq} -K2{iq}))/epsilon;
sdrbfKern.inverseWidth(indexIW) = invW;
%%% SWITCHING POINTS
indexSP = 3;
sptime = sdrbfKern.switchingTimes(indexSP);
sdrbfKern.switchingTimes(indexSP) = sptime + epsilon;
K1 = sdrbfKernCompute(sdrbfKern, t1, t1);
sdrbfKern.switchingTimes(indexSP) = sptime - epsilon;
K2 = sdrbfKernCompute(sdrbfKern, t1, t1);
sp1n = 0.5*sum(sum(K1{1} -K2{1}))/epsilon;
sdrbfKern.switchingTimes(indexSP) = sptime;
%%%
g = sdrbfKernGradient(sdrbfKern, t1, t1, covGrad);





