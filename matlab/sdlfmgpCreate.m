function  model = sdlfmgpCreate(q, d, X, y, options)

% SDLFMGPCREATE creates a switching dynamical latent force GP model. 
% FORMAT
% DESC Returns a structure for the switching dynamical latent force 
% Gaussian process model. This model corresponds to a multiple-output GP 
% in which each output is segmented in intervals and the intervals are 
% represented with localized latent force models. The continuity of the 
% outputs is kept forcing the outputs of the next interval to have initial 
% conditions given by the last point of the last interval.
% RETURN model : the structure for the sdlfmgp model
% ARG q : input dimension size.
% ARG d : output dimension size.
% ARG Y : set of training observations
% ARG X : set of training inputs
% ARG options : contains the options for the MULTIGP model, which includes
% if the model is approximated or full, the number of latent functions, the
% number of output functions.
%
% SEE ALSO : multigpCreate.m, gpCreate.m
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% SDLFMGP

if iscell(X)
    if size(X{end}, 2) ~= q
        error(['Input matrix X does not have dimension ' num2str(q)]);
    end
else
    if size(X, 2) ~= q
        error(['Input matrix X does not have dimension ' num2str(q)]);
    end
end
if iscell(y)
    % ensure it is a row vector of cells.
    y = y(:)';
    if size(y, 2) ~= d
        error(['Target cell array Y does not have dimension ' num2str(d)]);
    end
    for i = 1:size(y, 2)
        if(size(y{i}, 2)>1)
            error('Each element of the cell array should be a column vector.')
        end
    end
else
    if size(y, 2)~=d
        error(['Target matrix Y does not have dimension ' num2str(d)]);
    end
end


model.type = 'sdlfmgp';
model.q = q;

% Short hand for the kernel type
model.kernType = options.kernType;

% Number of latent functions
model.nlfPerInt = options.nlfPerInt;

% Number of intervals
model.nIntervals = options.nIntervals;

% Number of output functions
model.d = d;
model.nout = model.d;
model.numPositions = model.nout/(1+options.includeVel + options.includeAccel);
model.approx = options.approx;

% Set up default scale and bias for outputs
if isfield(options, 'scale') && ~isempty(options.scale)
    model.scale = options.scale;
else
    model.scale = ones(1, model.d);
end
if isfield(options, 'bias') && ~isempty(options.bias)
    model.bias = options.bias;
else
    model.bias = zeros(1, model.d);
end

% Initialization of the inputs for the model

switch model.approx
    case 'ftc'
        model.X = X;
        model.y = [];
        for i = 1:length(y)
            model.y = [model.y; y{i}];
        end
        model.N = size(model.y,1);
    case {'dtc','fitc','pitc', 'dtcvar'}
        % By default we assume all latent forces have the same number of
        % inducing points. For multiresolution models this might not be the
        % case.
        numActive = options.numActive;
        posX = zeros(numActive, q);
        switch options.initialInducingPositionMethod
            case 'espaced'
                % Especially useful for 1D case. It allows to move the
                % pseudo-inputs a further to the original input range
                for j = 1:q,
                    factor = 0.05;
                    med = (max(X{1}(:,j)) - min(X{1}(:,j)))/2;
%                     posX(:,j) = linspace(min(X{1}(:,j)) - factor*med, ...
%                         max(X{1}(:,j)) + factor*med, numActive)';
                    posX(:,j) = linspace(options.kern.switchingTimes(1)+0.1, ...
                        max(X{1}(:,j)) + factor*med, numActive)';
                end
            case 'espacedInRange'
                % Especially useful for 1D case. It restricts the
                % initial position of pseudo-inputs to be within the input range
                for j = 1:q,
                    posX(:,j) = linspace(min(X{1}(:,j)), max(X{1}(:,j)), numActive)';
                end
            case 'randomDataIsotopic'
                % The initial positions of pseudo-inputs are taken from
                % the data. In the isotopic case, since all inputs all
                % equal for each output, we use only X(1) to select the
                % positions, as long as the number of numActive is less
                % than the size of inputs.
                if size(X{1},1) >= numActive,
                    totX = cell2mat(X(1));
                else
                    totX = cell2mat(X');
                end
                index = randperm(size(totX,1));
                posX = totX(index(1:numActive),:);
            case 'randomDataHeterotopic'
                totX = cell2mat(X');
                index = randperm(size(totX,1));
                posX = totX(index(1:numActive),:);
            case 'randomComplete'
                posX = 0.5*rand(numActive, q);
            case 'fixIndices'
                posX = X{1}(options.fixIndices,:);
            case 'kmeansIsotopic'
                posX = kmeanlbg(X{1},numActive);
            case 'kmeansHeterotopic'
                totX = cell2mat(X');
                posX = kmeanlbg(totX,numActive);
            otherwise
                error('This is not valid initialization method for the input variables');
        end
        model.X{1,1} = posX;
        model.y = [];
        for i = 1:length(y)
            model.y = [model.y; y{i}];
            model.X{i + 1,1} = X{i};
            model.sizeX(i) = size(X{i},1);
        end
        model.N = size(model.y,1);
        X = model.X;
    otherwise
        error('Unknown model approximation')
end

model.includeInd = options.includeInd;
model.tieIndices = options.tieOptions.tieIndices;
model.includeNoise = options.includeNoise;
kernType = cell(1, 1 + options.includeNoise + options.includeInd);
cont = 0;
cont = cont + 1;
kernType{cont} = sdlfmgpKernComposer(options.kernType, model.nout, model.approx, options);

% To include independent kernel
if model.includeInd
    cont = cont + 1;
    kernType{cont} = multigpKernComposer('rbf', model.d, model.nlf, model.approx);
end
% To include noise
if model.includeNoise && strcmp(model.approx, 'ftc')
    cont = cont + 1;
    kernType{cont} = multigpKernComposer('white', model.d, 0, model.approx);
end

model.kern = kernCreate(X, {'cmpnd', kernType{:}});

% Additional kernel options or modifications
model.kern = sdmultiKernParamInit(model.kern, options);

switch model.approx
    case 'ftc'
        % Here we modify the multi kern that contains the parameters of the
        % switching dynamical kernel and create the structure for the
        % switching dynamical rbf kernel
           model.sdrbfKern = kernCreate(X, {'parametric', options.kern, 'sdrbf'});
        
    case {'fitc', 'pitc'}
        

end

if isfield(options, 'optimiser') && ~isempty(options.optimiser)
    model.optimiser = options.optimiser;
end

% NEIL: Learning of scales hasn't been included, although it should be.
model.learnScales = options.learnScales;
model.scaleTransform = optimiDefaultConstraint('positive');

switch model.approx
    case {'dtc','fitc','pitc', 'dtcvar'}
        model = spsdlfmgpCreate( model, options);
end

model.nParams = 0;

% Set up a mean function if one is given.
if isfield(options, 'meanFunction') && ~isempty(options.meanFunction)
    if isstruct(options.meanFunction)
        model.meanFunction = options.meanFunction;
    else
        if ~isempty(options.meanFunction)
            model.meanFunction = meanCreate(q, model.nout, X, y, options.meanFunctionOptions);
        end
    end
    model.nParams = model.nParams + model.meanFunction.nParams;
end

% Tie options according to the particular kernel employed
fhandle = [model.kernType 'SdlfmgpTieParam'];
if exist(fhandle, 'file')
    fhandle = str2func(fhandle);
    tieInd = fhandle(model, options);
else
    error('Function for tying parameters for this kernel type not implemented yet')
end

model.nParams = model.nParams + model.kern.nParams;

% Creates the noise model (Default model: one beta per output)
switch model.approx
    case {'dtc','fitc','pitc', 'dtcvar'}
        % In this structure is easier to put noise in the latent functions
        % at the top level
        if isfield(options, 'gamma') && ~isempty(options.gamma)
            if size(options.gamma,2) == model.nlfPerInt
                model.gamma = options.gamma;
            else
                model.gamma = options.gamma*ones(1,model.nlfPerInt);
            end
            model.gammaTransform =  optimiDefaultConstraint('positive');
            model.nParams = model.nParams + model.nlfPerInt;
        end
        if isfield(options, 'beta') && ~isempty(options.beta)
            if size(options.beta,2) == model.nout
                model.beta = options.beta;
            else
                model.beta = options.beta*ones(1,model.nout);
            end
            model.betaTransform =  optimiDefaultConstraint('positive');
            model.nParams = model.nParams + model.nout;      
        end
        if ~options.fixInducing
            model.nParams = model.nParams + sum(model.k)*model.q;
        end
end



fhandle = [model.kernType 'SdlfmgpFixParam'];
if exist(fhandle, 'file')
    fhandle = str2func(fhandle);
    model = fhandle(model);
end

model = modelTieParam(model, tieInd);
params = modelExtractParam(model);
model = modelExpandParam(model, params);

model.alpha = [];
end
