function [params, names] = spsdlfmgpExtractParam(model, paramPart, names)

% SPSDLFMGPEXTRACTPARAM Extract a parameter vector from a sparse SDLFMGP model.
% FORMAT
% DESC extracts the model parameters from a structure containing
% the information about a sparse swicthing dynamical LFM Gaussian process.
% ARG model : the model structure containing the information about
% the model.
% ARG paramPart: the parameters corresponding to the sparse part of the 
% multigp model
% RETURN params : a vector of parameters from the model.
%
% DESC extracts the model parameters from a structure containing
% the information about a sparse swicthing dynamical LFM Gaussian process.
% ARG model : the model structure containing the information about
% the model.
% ARG paramPart: the parameters corresponding to the sparse part of the 
% multigp model
% ARG names : cell array correspondig to the names of the basic multigp
% model structure.
% RETURN params : a vector of parameters from the model.
% RETURN names : cell array of parameter names.
%
% SEEALSO : multigpCreate, multigpExpandParam, modelExtractParam
%
% COPYRIGHT : Mauricio A Alvarez, 2008
%
% MODIFICATIONS : Mauricio A Alvarez, 2009, 2010

% SDLFMGP

switch model.approx
    case {'dtc','fitc', 'pitc', 'dtcvar'}
        if model.fixInducing
            params = paramPart;
            model.X_u = model.X{1};
            for k=2:model.nlf
                model.X_u = [model.X_u; model.X{k}];
            end
        else
            model.X_u = model.X{1};
            for k=2:model.nlf
                model.X_u = [model.X_u; model.X{k}];
            end
            params =  [paramPart model.X_u(:)'];            
            if nargout>1
                X_uNames = cell(size(model.X_u));
                if exist('Xunames.txt', 'file')
                    fidNames = fopen('Xunames.txt','r');
                    for i = 1:size(model.X_u, 1)
                        for j = 1:size(model.X_u, 2)
                            X_uNames{i, j} = fgetl(fidNames);
                        end
                    end
                    fclose(fidNames);
                else
                    for i = 1:size(model.X_u, 1)
                        for j = 1:size(model.X_u, 2)
                            X_uNames{i, j} = ['X_u (' num2str(i) ', ' num2str(j) ')'];
                        end
                    end
                end
                names = {names{:}, X_uNames{:}};
            end
        end
    otherwise
        %
end