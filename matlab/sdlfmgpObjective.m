function f = sdlfmgpObjective(params, model)

% SDLFMGPOBJECTIVE Wrapper function for SDLFMGPOPTIMISE objective.

% SDLFMGP

% We can use the same likelihood of the multigp model. 

model = modelExpandParam(model, params);
f = - sdlfmgpLogLikelihood(model);