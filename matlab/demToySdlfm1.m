% DEMTOYSDLFM Demonstrate switching dynamical latent force model on CMU data.

% SDLFMGP
clc
clear
randn('state', 1e6)
rand('twister', 1e6)

dataSetName = 'toySdlfm1';
experimentNo = 1;


[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);

yTemp2 = yTemp;
XTemp2 = XTemp;
yTemp = yTemp(1);
XTemp = XTemp(1);
% Just to test the code


options = sdlfmgpOptions('ftc');
options.nIntervals = 3;
options.nlfPerInt = 1;
options.kern.nIntervals = options.nIntervals;
options.kern.switchingTimes = [-1 7.6 6.6];
options.kern.nlfPerInt = options.nlfPerInt;
% options.optimiser = 'conjgrad';

Y = cell(1, length(yTemp));
X = cell(1, length(yTemp));

for i=1:length(yTemp)
%     options.bias(i) = mean(yTemp{i});
%    options.scale(i) = std(yTemp{i});
    Y{i} = yTemp{i};
    X{i} = XTemp{i};
end

% Get the time index.

q = 1;
d = length(yTemp);

% Creates the model

model = sdlfmgpCreate(q, d, X, Y, options);

% % Change values of initial parameters
params = modelExtractParam(model);
indexDamper = paramNameRegularExpressionLookup(model, '.* damper');
params(indexDamper) = log(0.4);
indexSpring = paramNameRegularExpressionLookup(model, '.* spring');
params(indexSpring) = log(1);
indexIW = paramNameRegularExpressionLookup(model, '.* inverse width .*');
params(indexIW) = log(1e-3);
indexVariance = paramNameRegularExpressionLookup(model, '.* variance');
params(indexVariance) = log(10);
model = modelExpandParam(model, params);

display = 1;
iters = 10;

% Trains the model and counts the training time
model = modelOptimise(model, [], [], display, iters);

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

% model.kern.comp{1}.comp{1}.sensitivity = (options.scale(1)^2)*model.kern.comp{1}.comp{1}.sensitivity;
% params = modelExtractParam(model);
% model = modelExpandParam(model, params);

% ll = 0;
% for j=1:length(yTemp2)
%     %yTest2 = zscore(yTemp2{j});
% %     yTest2 = yTemp2{j} - mean(yTemp2{j});
%     yTest2 = yTemp2{j};
%     ll = ll - model.N*log(2*pi) - model.logDetK - yTest2'*model.invK*yTest2;
% end

[XGT, void, void, fGT] = mapLoadData(dataSetName);
sdlfmgpToyResults(dataSetName, experimentNo, XTemp, yTemp, ...
    XGT, fGT)