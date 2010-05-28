% DEMTOYSDLFM Demonstrate switching dynamical latent force model on CMU data.

% SDLFM

clc
clear
randn('state', 1e6)
rand('twister', 1e6)

dataSetName = 'toySdlfm2';
experimentNo = 2;

[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);


indexSample = 1;
yTemp2 = yTemp;
XTemp2 = XTemp;
yTemp = yTemp{indexSample};
XTemp = XTemp{indexSample};

options = sdlfmgpOptions('dtcvar');
options.nIntervals = 2;
options.nlfPerInt = 1;
options.nlf = options.nlfPerInt;
options.kern.nIntervals = options.nIntervals;
options.kern.switchingTimes = [-1 5];
options.kern.nlfPerInt = options.nlfPerInt;
options.numActive = 20;
options.includeNoise = false;
options.beta = 1e-3;
options.initialInducingPositionMethod = 'espaced';


Y = cell(1, length(yTemp));
X = cell(1, length(yTemp));
scaleVal = zeros(1, length(yTemp));

for i=1:length(yTemp)
    Y{i} = yTemp{i};
    X{i} = XTemp{i};
%     options.bias(i)  = mean(Y{i});
    options.scale(i) = std(Y{i}); 
%     scaleVal(i) = var(Y{i});
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
% indexIW = paramNameRegularExpressionLookup(model, '.* inverse width .*');
% %params(indexIW) = log(1e-3*rand(1,length(indexIW)));
% params(indexIW) = log(1e-3);
% sens = paramNameRegularExpressionLookup(model, '.* sensitivity .*');
% params(sens) = [10 1 10 5 -10 1];
%params(sens) = 2*rand(1,length(sens));
model = modelExpandParam(model, params);

display = 1;
%if options.nIntervals == 1
iters = 100;
model = modelOptimise(model, [], [], display, iters);
% else
    % Trains the model and counts the training time
%     iters = 25;
%     model = fixDynSwitchingPointsParams(model, model2, false);
%     model = modelOptimise(model, [], [], display, iters);
%     fixN = length(model.fix);
%     model.fix(21:fixN) = [];
%     iters = 10;
%     model = fixDynSysParams(model, model2, false);
%     model = modelOptimise(model, [], [], display, iters);
%     fixN = length(model.fix);
%     model.fix(21:fixN) = [];
%     model = fixDynSwitchingPointsParams(model, model2, false);
%     model = modelOptimise(model, [], [], display, iters);
%     fixN = length(model.fix);
%     model.fix(21:fixN) = [];
%     model = fixDynSysParams(model, model2, false);
%     model = modelOptimise(model, [], [], display, iters);
%     fixN = length(model.fix);
%     model.fix(21:fixN) = [];
%     iters = 15;
%     model = modelOptimise(model, [], [], display, iters);
% end

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

sdlfmgpToyResults(dataSetName, experimentNo, XTemp2{indexSample}, yTemp2{indexSample}, ...
    XTestTemp{indexSample}, yTestTemp{indexSample})