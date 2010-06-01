function model = sdlfmvSdlfmgpFixParam(model)

% SDLFMVSDLFMGPFIXPARAM Fix parameters for a sdlfmgp gp with sdlfmv kernel
% FORMAT
% DESC Fix the parameters for a sdlfmgp model that uses SDLFMV kernel.
% RETURN model : model with fixed parameters included
% ARG model    : model before fixing the parameters
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% SDLFMGP

model = sdlfmSdlfmgpFixParam(model);