function  ll = sdlfmgpLogLikelihood( model)

% SDLFMGPLOGLIKELIHOOD Compute the log likelihood of a SDLFMGP.
% FORMAT
% DESC
% ARG model : sdlfmgp model structure for which the likelihood will be
% computed
% RETURN ll : value of the loglikelihood
%
% COPYRIGHT : Mauricio A. Alvarez, 2008
%
% MODIFICATIONS : Mauricio A. Alvarez, 2010 

% SDLFMGP

switch model.approx
    case 'ftc'
        dim = size(model.m, 1);
        ll = -dim*log(2*pi) - model.logDetK - model.m'*model.invK*model.m;
        ll = ll*0.5;
    case {'dtc','fitc','pitc','dtcvar'}        
        ll = spsdlfmgpLogLikelihood( model);
end
    
