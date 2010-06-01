function model = spsdlfmgpCreate(model, options)

% SPSDLFMGPCREATE
% DESC incorporates into de model options retaled with the sparse methods
% RETURN model : the structure for the sparse multigp model
% ARG model : input sparse model.
% ARG options : contains the options for the sparse multigp model

% COPYRIGHT : Mauricio Alvarez  2008

% MODIFICATIONS : Mauricio Alvarez 2009, 2010

% SDLFMGP

switch options.approx
    case {'dtc','fitc', 'pitc', 'dtcvar'}
        % Sub-sample inducing variables.
        model.k = options.numActive;
        model.fixInducing = options.fixInducing;
        model.X_u = model.X{1};
        for k=2:options.nlfPerInt
            model.X_u = [model.X_u; model.X{k}];
        end        
end
if sum(model.k)>model.N
    error('Number of active points cannot be greater than number of data.')
end
