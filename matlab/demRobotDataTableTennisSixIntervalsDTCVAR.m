% DEMCMU07WALKFEET Demonstrate sd latent force model on CMU data.

% SDLFMGP

clc
clear

rand('seed', 1e6);
randn('seed', 1e6);

addToolboxes(0,1)
% load data
P =8;
lengthSignal = 21;
load(['./robotData/robotDataP' num2str(P) 'Length' num2str(lengthSignal) '.mat'], 'X', 'y');


options = sdlfmgpOptions('dtcvar');
options.nIntervals = 6;
options.nlfPerInt = 1;
options.kern.nIntervals = options.nIntervals;
options.kern.switchingTimes = [-1 5 3 3 3 3];
options.kern.nlfPerInt = options.nlfPerInt;
options.kern.isNormalised = false;
if ~strcmp(options.approx,'ftc')
    options.numActive = 50;
    options.beta = 1e-1;
    options.gamma = 0;
    options.initialInducingPositionMethod = 'espaced';
    options.includeNoise = false;
else
    options.includeNoise = true;
end

y = y(1:4);
X = X(1:4);

for i=1:length(y)
    options.scale(i) = std(y{i});
    options.bias(i) = mean(y{i});
end


q = 1;
d = size(y,2);

% Creates the model

model = sdlfmgpCreate(q, d, X, y, options);
model.fix = [];

% We fix some additional parameters, namely, the value of the first
% switching point at -1 and the value of the covariances for the initial 
% conditions. 
model2 = model;
model2.nParams = size(model2.paramGroups,1);
model2.paramGroups = speye(model2.nParams);
if isfield(model, 'fix') && ~isempty(model.fix)
    count = length(model.fix);
else
    count = 0;
end
index = paramNameRegularExpressionLookup(model2, '.* switching point interval 1');
for k=1:length(index);
    count = count + 1;
    model.fix(count).index = index(k);
    model.fix(count).value = model.kern.comp{1}.comp{1}.switchingTimes(1);
end

if ~strcmp(options.approx,'ftc')
    ts = diff(model.X{2});
    ts = ts(1);
    % Reorganize positions seudo-inputs
    X_u = linspace(model.X{2}(1)-10*ts, ...
        model.X{2}(end)+10*ts, options.numActive)';
    model.X_u = X_u;
    model.X{1} = X_u;
end

% Change values of initial parameters
params = modelExtractParam(model);
sens = paramNameRegularExpressionLookup(model, '.* sensitivity .*');
params(sens) = 0.1*randn(1, length(sens));
model = modelExpandParam(model, params);


display = 1;
iters = 1000;
model = modelOptimise(model, [], [], display, iters);
save('demRobotSixIntervalDTCVAR.mat', 'model');
