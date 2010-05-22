% DEMTOYSDLFM Demonstrate switching dynamical latent force model on CMU data.

% SDLFM

randn('state', 1e6)
rand('twister', 1e6)

%dataSetName = 'cmu49BalanceArm';
clear
dataSetName = 'toySdlfm1';
experimentNo = 1;


% load data
%[y, void, yTest, void] = lvmLoadData(dataSetName);
[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);

% Just to test the code
indexShuffled = [4 5 1 2 3]; 


yTemp = yTemp(indexShuffled);

yP = cell2mat((yTemp(1:4))');
sizeX = cellfun('length', yTemp(1:4));
scaleVal = std(yP);
yTemp2 = mat2cell(yP/scaleVal, sizeX, 1)';

options.type = 'sdlfmgp';
options.numModels = 4;
options.compOptions = sdlfmgpOptions('ftc');
options.compOptions.nIntervals = 3;
options.compOptions.nlfPerInt = 1;
options.compOptions.kern.nIntervals = options.compOptions.nIntervals;
options.compOptions.kern.switchingTimes = [-1 7.6 6.6];
options.compOptions.kern.nlfPerInt = options.compOptions.nlfPerInt;
options.separate = [];
options.optimiser = 'scg';

X{1}{1} = XTemp{1};  
X{2}{1} = XTemp{1};  
X{3}{1} = XTemp{1};  
X{4}{1} = XTemp{1};  
Y{1}{1} = yTemp2{1};
Y{2}{1} = yTemp2{2};
Y{3}{1} = yTemp2{3};
Y{4}{1} = yTemp2{4};


q = 1;
d = 1;

% Creates the model
model = multimodelCreate(q, d, {X{1}; X{2}; X{3}; X{4}}, {Y{1}; Y{2}; Y{3}; Y{4}}, options);

%model = sdlfmgpCreate(q, d, X, Y, options);

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
iters = 5;

% Trains the model and counts the training time
model = modelOptimise(model, [], [], display, iters);

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');
model2 = model;
model = model2.comp{1};
ll = 0;
for j=1:length(yTemp2)
    %yTest2 = zscore(yTemp2{j});
%     yTest2 = yTemp2{j} - mean(yTemp2{j});
    yTest2 = yTemp2{j};
    ll = ll - model.N*log(2*pi) - model.logDetK - yTest2'*model.invK*yTest2;
end
