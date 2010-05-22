% DEMTOYSDLFM Demonstrate switching dynamical latent force model on CMU data.

% SDLFM

randn('state', 1e6)
rand('twister', 1e6)

%dataSetName = 'cmu49BalanceArm';
dataSetName = 'toySdlfm1';
experimentNo = 2;

[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);

yTemp2 = yTemp;
XTemp2 = XTemp;
yTemp = yTemp(1);
XTemp = XTemp(1);


% Just to test the code

% Get the time index.
%fps = 120/32;
%scaleVal = sqrt(sum(var(y)));
%y = y/scaleVal;

options = sdlfmgpOptions('pitc');
options.nIntervals = 2;
options.nlfPerInt = 1;
options.nlf = options.nlfPerInt;
% options.includeVel = true;
% options.kern.includeVel  = true;
options.kern.nIntervals = options.nIntervals;
options.kern.switchingTimes = [-2  5];
options.kern.nlfPerInt = options.nlfPerInt;
options.numActive = 20;
options.includeNoise = false;
options.beta = 1e-3;
%options.meanFunction = true;
%options.meanFunctionOptions.type = 'lfm';
%options.meanFunctionOptions.nlfPerInt  = options.kern.nlfPerInt;
%options.meanFunctionOptions.nIntervals = options.kern.nIntervals;
options.optimiser = 'conjgrad';
options.initialInducingPositionMethod = 'espaced';

Y = cell(1, length(yTemp));
X = cell(1, length(yTemp));

for i=1:length(yTemp)
%     options.bias(i) = yTemp{i}(1);
    options.scale(i) = std(yTemp{i});
    Y{i} = yTemp{i};
    X{i} = XTemp{i};
end

%Y = cell(1, size(y,2));
%X = cell(1, size(y,2));

% for i = 1:size(y, 2)
%   Y{i} = y(1:35, i);
%   X{i} = (1:35)'/fps;
% end

% Get the time index.

q = 1;
d = length(yTemp);

% Creates the model

model = sdlfmgpCreate(q, d, X, Y, options);

%model = modelExpandParam(model, param);

display = 1;
iters = 50;

% Trains the model and counts the training time
model = modelOptimise(model, [], [], display, iters);

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

[XGT, void, void, fGT] = mapLoadData(dataSetName);
sdlfmgpToyResults(dataSetName, experimentNo, XTemp, yTemp, ...
    XGT, fGT)