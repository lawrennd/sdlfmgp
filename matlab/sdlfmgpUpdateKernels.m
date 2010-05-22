function model = sdlfmgpUpdateKernels(model)

% SDLFMGPUPDATEKERNELS Updates the kernel in a SDLFMGP structure.
% FORMAT
% DESC computes the kernel matrix for the switching dynamical latent force
% gp model. 
% ARG model : model containing the parameters used to compute the kernel
% RETURN model : model containing the kernel matrix and its inverse
%
% SEE ALSO : multigpUpdateKernels
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% SDLFMGP

switch model.approx
    case 'ftc'
        K = real(kernCompute(model.kern, model.X));
        model.K = K;
        [model.invK, U, jitter] = pdinv(K);
        if any(jitter>1e-4)
            fprintf('Warning: UpdateKernels added jitter of %2.4f\n', jitter)
        end
        model.alpha = multigpComputeAlpha(model);
        model.logDetK = logdet(model.K, U);
    case {'dtc','fitc','pitc', 'dtcvar'}
        model = spsdlfmgpUpdateKernels(model);
end


