% DEMTOYSDLFM1 Demonstrate switching dynamical latent force model on
% artificial data 

% SDLFMGP

clc
clear
randn('state', 1e6)
rand('twister', 1e6)

dataSetName = 'toySdlfm3';
experimentNo = 2;

[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);


indexSample = 1;
yTemp2 = yTemp;
XTemp2 = XTemp;
yTemp = yTemp{indexSample};
XTemp = XTemp{indexSample};


options = sdlfmgpOptions('dtcvar');
options.nIntervals = 1;
options.nlfPerInt = 2;
options.nlf = options.nlfPerInt;
options.kern.nIntervals = options.nIntervals;
options.kern.switchingTimes = [-2 10];
options.kern.nlfPerInt = options.nlfPerInt;
options.kern.isNormalised = true;
options.numActive = 50;
options.includeNoise = false;
options.beta = 1e-1;
options.gamma = exp(-2);
options.initialInducingPositionMethod = 'espaced';


Y = cell(1, length(yTemp));
X = cell(1, length(yTemp));
scaleVal = zeros(1, length(yTemp));

for i=1:length(yTemp)
    Y{i} = yTemp{i};
    X{i} = XTemp{i};
  %  options.bias(i)  = mean(Y{i});
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

% Reorganize positions seudo-inputs
% X_u = linspace(options.kern.switchingTimes(1)+2, ...
%     model.X{2}(end)+0.5, options.numActive)';
% 
% model.X_u = X_u;
% model.X{1} = X_u;

% Change values of initial parameters
params = modelExtractParam(model);
indexDamper = paramNameRegularExpressionLookup(model, '.* damper');
params(indexDamper) = log([0.4 1 1]);
indexSpring = paramNameRegularExpressionLookup(model, '.* spring');
params(indexSpring) = log([2 3 0.5]);
% indexIW = paramNameRegularExpressionLookup(model, '.* inverse width .*');
% params(indexIW) = log([1e-3 1]);
% % sens = paramNameRegularExpressionLookup(model, '.* sensitivity .*');
% % temp  = repmat(options.scale',1,options.nIntervals)';
% % scaleVector = temp(:)';
% % params(sens) = [1 5 5 1 1 1]./scaleVector;
% % sens = paramNameRegularExpressionLookup(model, '.* sensitivity .*');
% % params(sens) = [1 5];
% % indexVariance = paramNameRegularExpressionLookup(model, '.* variance');
% % params(indexVariance) = log(10);
model = modelExpandParam(model, params);


display = 1;
% % if options.nIntervals == 1
% %baseFix = length(model.fix);
% %model = fixDynSwitchingPointsParams(model, model2, false);
iters = 100;
model = modelOptimise(model, [], [], display, iters);
% else
%     iters = 20;
%     baseFix = length(model.fix);
%     model = fixDynSwitchingPointsParams(model, model2, false);
%     model = modelOptimise(model, [], [], display, iters);
%     fixN = length(model.fix);
%     model.fix(baseFix+1:fixN) = [];    
%     iters = 5;
%     numberRepeats =20 ;
%     for j=1:numberRepeats
%         model = fixDynSwitchingPointsParams(model, model2, false);
%         model = modelOptimise(model, [], [], display, iters);
%         fixN = length(model.fix);
%         model.fix(baseFix+1:fixN) = [];        
%         model = fixDynSysParams(model, model2, false);
%         model = modelOptimise(model, [], [], display, iters);
%         fixN = length(model.fix);
%         model.fix(baseFix+1:fixN) = [];        
%     end
%     iters = 10;
%     model = modelOptimise(model, [], [], display, iters);
% end

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

sdlfmgpToyResults(dataSetName, experimentNo, XTemp2{indexSample}, yTemp2{indexSample}, ...
    XTestTemp{indexSample}, yTestTemp{indexSample})