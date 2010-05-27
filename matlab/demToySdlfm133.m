% DEMTOYSDLFM133 Demonstrate switching dynamical latent force model on
% artificial data 

% SDLFMGP

clc
clear
randn('state', 1e6)
rand('twister', 1e6)
addToolboxes(0,1)

dataSetName = 'toySdlfm2';
experimentNo = 1;

[XTemp2, yTemp2, XTestTemp, yTestTemp] = mapLoadData(dataSetName);


options = sdlfmgpOptions('ftc');
options.nIntervals = 3;
options.nlfPerInt = 1;
options.kern.nIntervals = options.nIntervals;
options.kern.switchingTimes = [-1 5 8];
options.kern.nlfPerInt = options.nlfPerInt;


for indexSample = 1:length(yTemp2)
fprintf('TEST NUMBER : %d\n', indexSample);
yTemp = yTemp2{indexSample};
XTemp = XTemp2{indexSample};
Y = cell(1, length(yTemp));
X = cell(1, length(yTemp));
scaleVal = zeros(1, length(yTemp));
for i=1:length(yTemp)
    Y{i} = yTemp{i};
    X{i} = XTemp{i};
    options.bias(i)  = mean(Y{i});
end

% Get the time index.

q = 1;
d = length(yTemp);

% Creates the model

model = sdlfmgpCreate(q, d, X, Y, options);

% We fix some additional parameters, namely, the value of the first
% switching point at -1 and the value of the covariances for the initial 
% conditions. 
model2 = model;
model2.nParams = size(model2.paramGroups,1);
model2.paramGroups = speye(model2.nParams);
count = length(model.fix);
% Fix the value of the initial condition first
% Now, fix the covariances of the initial conditions
valkInit = model.kern.comp{1}.LIC(:)';
index = paramNameRegularExpressionLookup(model2, '.* kInit.*');
for k=1:length(index);
    count = count + 1;
    model.fix(count).index = index(k);
    model.fix(count).value = valkInit(k);
end
index = paramNameRegularExpressionLookup(model2, '.* switching point interval 1');
for k=1:length(index);
    count = count + 1;
    model.fix(count).index = index(k);
    model.fix(count).value = model.kern.comp{1}.comp{1}.switchingTimes(1);
end
% Change values of initial parameters
params = modelExtractParam(model);
indexDamper = paramNameRegularExpressionLookup(model, '.* damper');
params(indexDamper) = log([0.4 1]);
indexSpring = paramNameRegularExpressionLookup(model, '.* spring');
params(indexSpring) = log([2 3]);
%indexIW = paramNameRegularExpressionLookup(model, '.* inverse width .*');
%params(indexIW) = log(1e-2);
%sens = paramNameRegularExpressionLookup(model, '.* sensitivity .*');
%params(sens) = 1 + 5*rand(1,length(sens));
indexVariance = paramNameRegularExpressionLookup(model, '.* variance');
params(indexVariance) = log(20);
model = modelExpandParam(model, params);

display = 1;
iters = 100;
model = modelOptimise(model, [], [], display, iters);
modelBatch{indexSample} = model;
save(['toy1Batch' num2str(options.nIntervals) '.mat'], 'modelBatch');
end

