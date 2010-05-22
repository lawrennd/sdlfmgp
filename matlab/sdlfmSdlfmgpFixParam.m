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

%/~

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Just use temporarily to check something
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


index = paramNameRegularExpressionLookup(model, '.* switching point interval 1');

for k=1:length(index);
    count = count + 1;
    model.fix(count).index = index(k);
    model.fix(count).value = model.kern.comp{1}.comp{1}.switchingTimes(1);
end

valkInit = model.kern.comp{1}.LIC(:)';
index = paramNameRegularExpressionLookup(model, '.* kInit.*');

for k=1:length(index);
    count = count + 1;
    model.fix(count).index = index(k);
    model.fix(count).value = valkInit(k);
end
%~/
 