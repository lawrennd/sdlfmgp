function g = sdlfmgpGradient(params, model)

% SDLFMGPGRADIENT Gradient wrapper for a SDLFMGP model.

% SDLFMGP

model = modelExpandParam(model, params);
g = - modelLogLikeGradients(model);
