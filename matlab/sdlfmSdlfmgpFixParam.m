function model = sdlfmSdlfmgpFixParam(model)

% SDLFMSDLFMGPFIXPARAM Fix parameters for a sdlfmgp gp with sdlfm kernel
% FORMAT
% DESC Fix the parameters for a sdlfmgp model that uses SDLFM kernel.
% RETURN model : model with fixed parameters included
% ARG model    : model before fixing the parameters
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% SDLFMGP

% Fix the masses of the sdlfm kernels.
index = paramNameRegularExpressionLookup(model, '.* mass');
count = 0;
for k=1:length(index);
    count = count + 1;
    model.fix(count).index = index(k);
    model.fix(count).value = expTransform(0.1, 'xtoa');
end
 