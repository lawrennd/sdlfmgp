% DEMTOYSDLFM11 Demonstrate switching dynamical latent force model on
% artificial data

% SDLFMGP

randn('state', 1e6)
rand('twister', 1e6)

clear
dataSetName = 'toySdlfm1';
experimentNo = 1;

% load data
[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);

% Just to test the code
indexShuffled = randperm(10);
indexShuffled = randperm(10);
indexShuffled = randperm(10);

%indexShuffled = [1 2 9 3 4 5 6 7 8 10];

% yTemp = yTestTemp;
%yTemp = yTemp(indexShuffled(1:5))

nTotal = length(yTemp);
nTraining = 5;
nTesting = nTotal - nTraining;

yTemp = yTemp(indexShuffled);

% yP = cell2mat((yTemp(1:nTraining))');
% sizeX = cellfun('length', yTemp(1:nTraining));
% scaleVal = std(yP);
% yTemp2 = mat2cell(yP/scaleVal, sizeX, 1)';

yTemp2 = yTemp;

options.type = 'sdlfmgp';
options.numModels = nTraining;
options.separate = [];
options.optimiser = 'scg';

X = cell(1, nTraining);
Y = cell(1, nTraining);

for k=1:nTraining
   X{k}{1} = XTemp{1};
   Y{k}{1} = yTemp2{k};
   options.compOptions{k} = sdlfmgpOptions('ftc');
   options.compOptions{k}.nIntervals = 3;
   options.compOptions{k}.nlfPerInt = 1;
   options.compOptions{k}.kern.nIntervals = options.compOptions{k}.nIntervals;
   options.compOptions{k}.kern.switchingTimes = [-1 8 6];
   options.compOptions{k}.kern.nlfPerInt = options.compOptions{k}.nlfPerInt;
   options.compOptions{k}.kern.isNormalised = true;
%    options.compOptions{k}.bias = mean(yTemp2{k});
%    options.compOptions{k}.scale = std(yTemp2{k});
end

XT = { X{:} };
YT = { Y{:} };

q = 1;
d = 1;

% Creates the model
model = multimodelCreate(q, d, XT, YT, options);

% We fix some additional parameters, namely, the value of the first
% switching point at -1 and the value of the covariances for the initial 
% conditions. 
model2 = model.comp{1};
model2.nParams = size(model2.paramGroups,1);
model2.paramGroups = speye(model2.nParams);
% Fix the value of the initial condition first
for j=1:nTraining
    count = length(model.comp{j}.fix);
    % Fix the covariances of the initial conditions
    valkInit = model.comp{1}.kern.comp{1}.LIC(:)';
    index = paramNameRegularExpressionLookup(model2, '.* kInit.*');
    for k=1:length(index);
        count = count + 1;
        model.comp{j}.fix(count).index = index(k);
        model.comp{j}.fix(count).value = valkInit(k);
    end
    index = paramNameRegularExpressionLookup(model2, '.* switching point interval 1');
    for k=1:length(index);
        count = count + 1;
        model.comp{j}.fix(count).index = index(k);
        model.comp{j}.fix(count).value = model.comp{1}.kern.comp{1}.comp{1}.switchingTimes(1);
    end
end
% Change values of initial parameters
params = modelExtractParam(model);
indexDamper = paramNameRegularExpressionLookup(model, '.* damper');
params(indexDamper) = log(0.4);
indexSpring = paramNameRegularExpressionLookup(model, '.* spring');
params(indexSpring) = log(5);
% indexIW = paramNameRegularExpressionLookup(model, '.* inverse width .*');
% params(indexIW) = log(1e-5);
indexVariance = paramNameRegularExpressionLookup(model, '.* variance');
params(indexVariance) = log(10);
model = modelExpandParam(model, params);

display = 1;
tic
%if options.compOptions{1}.nIntervals == 1
iters = 50;
model = modelOptimise(model, [], [], display, iters);
%else
%     iters = 10;    
%     % Fix switching points only
%     model = fixDynSwitchingPointsParams(model, model2);
%     model = modelOptimise(model, [], [], display, iters);
%     iters = 10;
%     model = deleteFix(model, 7:8);
%     model = fixDynSysParams(model, model2);
%     model = modelOptimise(model, [], [], display, iters);
%     iters = 10;
%     model = deleteFix(model, 7:14);
%     model = fixDynSwitchingPointsParams(model, model2);
%     model = modelOptimise(model, [], [], display, iters);
%     model = deleteFix(model, 7:8);
%     model = fixDynSysParams(model, model2);
%     model = modelOptimise(model, [], [], display, iters);
%     model = deleteFix(model, 7:14);
%     model = modelOptimise(model, [], [], display, iters);
% end
toc
% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

% Use the same scaleVal from the training set to normalize the test set
% yP = cell2mat(yTemp(nTraining+1:nTotal)');
% sizeX = cellfun('length', yTemp(nTraining+1:nTotal));
% yTest2c = mat2cell(yP/scaleVal, sizeX, 1)';

yTest2c = yTemp(nTraining+1:nTotal);
ll = 0;
error = zeros(1, nTesting);
for j=1:nTesting
    yTest2 = yTest2c{j};
    [mu, varsigma] = sdlfmgpPosteriorMeanVar(model.comp{1}, XTemp{1});
    error(j) = mean((yTest2 - mu{2}).^2)/var(yTest2);
%     ll = ll - model.comp{1}.N*log(2*pi) - model.comp{1}.logDetK ...
%         - yTest2'*model.comp{1}.invK*yTest2;
end
% avgll = ll/nTesting;