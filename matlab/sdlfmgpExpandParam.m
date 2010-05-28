function model = sdlfmgpExpandParam(model, params)

% SDLFMGPEXPANDPARAM Expand the given parameters into a SDLFMGP structure.
% FORMAT
% DESC expands the model parameters from a structure containing
% the information about a switching dynamical latent force GP model.
% ARG model : the model structure containing the information about
% the GP model.
% ARG params : a vector of parameters from the model.
% RETURN model : the model structure containing the information about
% the model updated with the new parameter vector.
%
% SEEALSO : sdlfmgpCreate, modelExpandParam, modelExtractParam
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% SDLFMGP

paramPart = real(params);

if isfield(model, 'fix')
    for i = 1:length(model.fix)
       paramPart(model.fix(i).index) = model.fix(i).value;
    end
end


startVal = 1;
endVal = model.kern.nParams;
kernParams = paramPart(startVal:endVal);

if length(kernParams) ~= model.kern.nParams
    error('Kernel Parameter vector is incorrect length');
end

model.kern = kernExpandParam(model.kern, kernParams);

% Check if there is a mean function.
if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
    startVal = endVal + 1;
    endVal = endVal + model.meanFunction.nParams;
    model.meanFunction = meanExpandParam(model.meanFunction, ...
        paramPart(startVal:endVal));
end

% Check if there is a gamma parameter.
if isfield(model, 'gamma') && ~isempty(model.gamma)
    startVal = endVal + 1;
    endVal = endVal + model.nlf;
    fhandle = str2func([model.gammaTransform 'Transform']);
    model.gamma = fhandle(paramPart(startVal:endVal), 'atox');
end

% Check if there is a beta parameter.
if isfield(model, 'beta') && ~isempty(model.beta)
    startVal = endVal + 1;
    endVal = endVal + model.nout;
    fhandle = str2func([model.betaTransform 'Transform']);
    model.beta = fhandle(paramPart(startVal:endVal), 'atox');
end

switch model.approx
    case {'dtc','fitc','pitc', 'dtcvar'}
        if ~model.fixInducing
            startVal = endVal + 1;
            endVal = endVal +  sum(model.k)*model.q;
            model = spmultigpExpandParam(model, paramPart(startVal:endVal)); % must do the contrary of spmultigpExtractParam
        end
    otherwise
        %
end
% Keep the values of parameters related with the mean at the top level.
model.m = sdlfmgpComputeM(model);
model = sdlfmgpUpdateKernels(model);

% Update the vector 'alpha' for computing posterior mean.
if isfield(model, 'alpha')
    model.alpha = multigpComputeAlpha(model);
end
end