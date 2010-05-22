function model = spsdlfmgpUpdateKernels(model)

% SPSDLFMGPUPDATEKERNELS 
% FORMAT
% DESC Update the kernels that are needed for the sparse sdlfmgp
% ARG model : the model structure containing the model parameters
% RETURN model : the model structure with updated kernels.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% SDLFMGP

switch model.approx
    case {'dtc','fitc', 'pitc','dtcvar'}      
        model = spsdlfmgpKernCompute(model);    
    otherwise
        %
end

if isfield(model, 'beta') && ~isempty(model.beta)
    model = spmultigpUpdateAD(model);
end