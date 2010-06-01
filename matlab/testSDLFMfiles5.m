% TEST ALL THE DERIVATIVES FOR POSPOS, VELPOS, POSVEL, VELVEL, ACCELVEL,
% VELACCEL, ACCELACCEL

%clc
clear
%randn('state', 1e6)
%rand('twister', 1e6)
% numSeed = rand('twister');
% numSeed = numSeed(1);

randn('state', 1967688546)
rand('twister', 1967688546)

% Create sdlfmKern1
nIntervals = 2;
tInit = 0;
tFinal = 20;
nPoints = 200;

t1 = linspace(tInit,tFinal, nPoints)';


options.nIntervals = nIntervals;
options.switchingTimes = [-2  10];
options.nlfPerInt = 1;
options.isNormalised = true;

sdlfmKern1 = kernCreate(t1, {'parametric', options ,'sdlfm'});
[params, names] = kernExtractParam(sdlfmKern1);
inverseWidth = 2*rand(options.nlfPerInt, options.nIntervals);
sensitivity = 0.5 + rand(options.nlfPerInt, options.nIntervals);
% sdlfmKern1.inverseWidth = [1e-3 1];
sdlfmKern1.sensitivity = sensitivity;
sdlfmKern1.mass = 0.1;
sdlfmKern1.damper = 0.4;
sdlfmKern1.spring = 2;
% sdlfmKern1.sensitivity = [1 5];
% Create sdlfmKern2

sdlfmKern2 = kernCreate(t1, {'parametric', options ,'sdlfm'});
params = kernExtractParam(sdlfmKern2);
params(sdlfmKern2.outputIndx) = 2*rand(1,length(sdlfmKern2.outputIndx)); % Change the mass, damper, spring
sdlfmKern2.inverseWidth = inverseWidth;
sdlfmKern2 = kernExpandParam(sdlfmKern2, params);


% Create sdrbfKern

sdrbfKern = kernCreate(t1, {'parametric', options, 'sdrbf'});
sdrbfKern.inverseWidth = sdlfmKern1.inverseWidth;

K = sdlfmXsdrbfKernCompute(sdlfmKern1, sdrbfKern, t1, t1);


fhandle1 = str2func('sdlfmXsdrbfKernCompute');
fhandle2 = str2func('sdlfmXsdrbfKernGradient');

covGrad = ones(length(t1));
iq = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MASS
epsilon = 1e-6;
mass = sdlfmKern1.mass;
sdlfmKern1.mass = mass + epsilon;
K1 = fhandle1(sdlfmKern1, sdrbfKern, t1, t1);
sdlfmKern1.mass = mass - epsilon;
K2 = fhandle1(sdlfmKern1, sdrbfKern, t1, t1);
mass1n = 0;
for iq=1:options.nlfPerInt
    mass1n = mass1n + 0.5*sum(sum((K1{iq} - K2{iq})))/epsilon;
end
iq = 1;
%%%%%%%%%%%%%%%%%%%%%
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
        gradNum(i,j) = 0.5*(sum(sum(K1{iq}(start1:end1, start2:end2) - K2{iq}(start1:end1, start2:end2))))/epsilon;
        start2 = end2 + 1;
    end
    start1 = end1 + 1;
end

%%%%%%%%%%%%%%%%%%%%%

sdlfmKern1.mass = mass;
%%%% SPRING 1
spring = sdlfmKern1.spring;
sdlfmKern1.spring = spring + epsilon;
K1 = fhandle1(sdlfmKern1, sdrbfKern, t1, t1);
sdlfmKern1.spring = spring - epsilon;
K2 = fhandle1(sdlfmKern1, sdrbfKern, t1, t1);
spring1n = 0;
for iq=1:options.nlfPerInt
    spring1n = spring1n + 0.5*sum(sum((K1{iq} - K2{iq})))/epsilon;
end
sdlfmKern1.spring = spring;
%%%% DAMPER 1
damper = sdlfmKern1.damper;
sdlfmKern1.damper = damper + epsilon;
K1 = fhandle1(sdlfmKern1, sdrbfKern, t1, t1);
sdlfmKern1.damper = damper - epsilon;
K2 = fhandle1(sdlfmKern1, sdrbfKern, t1, t1);
damper1n = 0;
for iq=1:options.nlfPerInt
    damper1n = damper1n + 0.5*sum(sum((K1{iq} - K2{iq})))/epsilon;
end
sdlfmKern1.damper = damper;
%%%% INVERSE WIDTH 1
iq = 1;
indexIW = 1;
invW = sdlfmKern1.inverseWidth(indexIW);
sdlfmKern1.inverseWidth(indexIW) = invW + epsilon;
sdrbfKern.inverseWidth(indexIW) = invW + epsilon;
K1 = fhandle1(sdlfmKern1, sdrbfKern, t1, t1);
sdlfmKern1.inverseWidth(indexIW) = invW - epsilon;
sdrbfKern.inverseWidth(indexIW) = invW - epsilon;
K2 = fhandle1(sdlfmKern1, sdrbfKern, t1, t1);
inv1n = 0.5*sum(sum((K1{iq} - K2{iq})))/epsilon;
sdlfmKern1.inverseWidth(indexIW) = invW;
sdrbfKern.inverseWidth(indexIW) = invW;
%%%% SENSITIVITY 1
indexS = indexIW;
senS = sdlfmKern1.sensitivity(indexS);
sdlfmKern1.sensitivity(indexS) =  senS + epsilon;
K1 = fhandle1(sdlfmKern1, sdrbfKern, t1, t1);
sdlfmKern1.sensitivity(indexS) =  senS - epsilon;
K2 = fhandle1(sdlfmKern1, sdrbfKern, t1, t1);
sensn = 0.5*(sum(sum(K1{iq} - K2{iq})))/epsilon;
sdlfmKern1.sensitivity(indexS) =  senS;
%%%% SWITCHING TIME
indexSP = 1;
spTimes = sdlfmKern1.switchingTimes(indexSP);
sdlfmKern1.switchingTimes(indexSP) = spTimes + epsilon;
sdrbfKern.switchingTimes(indexSP) = spTimes + epsilon;
K1 = fhandle1(sdlfmKern1, sdrbfKern, t1, t1);
sdlfmKern1.switchingTimes(indexSP) = spTimes - epsilon;
sdrbfKern.switchingTimes(indexSP) = spTimes - epsilon;
K2 = fhandle1(sdlfmKern1, sdrbfKern, t1, t1);
sp1n = 0;
for iq=1:options.nlfPerInt
    sp1n = sp1n + 0.5*sum(sum((K1{iq} - K2{iq})))/epsilon;
end
sdlfmKern1.switchingTimes(indexSP) = spTimes;
sdrbfKern.switchingTimes(indexSP) = spTimes;


[g1, g2] = fhandle2(sdlfmKern1, sdrbfKern, t1, t1, covGrad);
stop = 1;






