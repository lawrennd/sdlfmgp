function tieInd = sdlfmvSdlfmgpTieParam(model, options)

% SDLFMVSDLFMGPTIEPARAM Tie parameters for a sdlfmgp model with sdlfmv kernel
% FORMAT
% DESC Tie the parameters for a sdlfmgp model that uses a sdlfm kernel.
% RETURN tieInd : cell with elements containing the indexes of parameters
% to tie.
% ARG model : model created
% ARG options : options for the model.
%
% COPYRIGHT : Mauricio A. Alvarez 2010

% SDLFMGP

tieInd = sdlfmSdlfmgpTieParam(model, options);




