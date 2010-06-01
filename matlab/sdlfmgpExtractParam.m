function [param, names] = sdlfmgpExtractParam(model)

% SDLFMGPEXTRACTPARAM Extract the parameters of a SDLFMGP model.
% FORMAT
% DESC extracts the model parameters from a structure containing
% the information about a switching dynamical latent force Gaussian process.
% ARG model : the model structure containing the information about
% the model.
% RETURN params : a vector of parameters from the model.
%
% DESC does the same as above, but also returns parameter names.
% ARG model : the model structure containing the information about
% the model.
% RETURN params : a vector of parameters from the model.
% RETURN names : cell array of parameter names.
%
% SEEALSO : sdlfmgpCreate, multigpExtractParam, modelExtractParam
%
% COPYRIGHT : Mauricio A Alvarez, 2010

% SDLFMGP

if nargout>1
    [kernParams, kernParamNames] = kernExtractParam(model.kern);  
else
    kernParams = kernExtractParam(model.kern);
end

kernParams = real(kernParams);

% Check if the output scales are being learnt.
if model.learnScales
    fhandle = str2func([model.scaleTransform 'Transform']);
    scaleParams = fhandle(model.scale, 'xtoa');
    scaleParamNames= cell(1,length(scaleParams));
    if nargout > 1
        for i = 1:length(scaleParams)
            scaleParamNames{i} = ['Output Scale ' num2str(i)];
        end
    end
else
    scaleParams = [];
    scaleParamNames = {};
end

% Check if there is a mean function.
if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
    if nargout>1
        [meanFuncParams, meanFuncParamNames] = meanExtractParam(model.meanFunction);
        for i = 1:length(meanFuncParamNames)
            meanFuncParamNames{i} = ['mean Func ' meanFuncParamNames{i}];
        end
    else
        meanFuncParams = meanExtractParam(model.meanFunction);
    end
else
    meanFuncParamNames = {};
    meanFuncParams =[];
end

if isfield(model, 'gamma') && ~isempty(model.gamma)
    fhandle = str2func([model.gammaTransform 'Transform']);
    gammaParams = fhandle(model.gamma, 'xtoa');    
    if nargout>1
        gammaParamNames = cell(model.nlfPerInt,1);
        for i = 1:length(gammaParams)
            gammaParamNames{i} = ['Gamma ' num2str(i)];
        end
    end
else
    gammaParamNames = {};
    gammaParams =[];
end

if isfield(model, 'beta') && ~isempty(model.beta)
    fhandle = str2func([model.betaTransform 'Transform']);
    betaParams = fhandle(model.beta, 'xtoa');    
    if nargout>1
        betaParamNames = cell(model.nout,1);
        for i = 1:length(betaParams)
            betaParamNames{i} = ['Beta ' num2str(i)];
        end
    end
else
    betaParamNames = {};
    betaParams =[];
end


% Check if there is a noise model
paramPart = [kernParams scaleParams meanFuncParams gammaParams betaParams];

if nargout > 1
    names = {kernParamNames{:}, meanFuncParamNames{:}, ...
        scaleParamNames{:}, gammaParamNames{:}, betaParamNames{:}};
end

if isfield(model, 'fix')
    for i = 1:length(model.fix)
        paramPart(model.fix(i).index) = model.fix(i).value;
    end
end


switch model.approx
    case 'ftc'
        param = paramPart;
    case {'dtc','fitc','pitc', 'dtcvar'}
        if nargout>1
            [param, names] = spmultigpExtractParam(model, paramPart, names);
        else
            param = spmultigpExtractParam(model, paramPart);
        end
    otherwise
        error('Unknown approximation')
end
