function model = spsdlfmgpKernCompute(model)

% SPSDLFMGPKERNCOMPUTE Computes the kernels for the sparse approximation.
% FORMAT  
% DESC computes the kernels for the sparse approximation in the SDLFM GP
% model. We decide how to store the covariances of the
% outputs: as a vector, corresponding to the diagonal of the
% covariances of the outputs for the fitc or as a full covariance,
% corresponding to the pitc.
% ARG model : the model structure 
% RETURN model : the modified model structure with the kernels updated.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% SDLFMGP


% Get the initial conditions

kyyTemp = model.kern.comp{1}.KIC(1:model.kern.comp{1}.numPositions, 1:model.kern.comp{1}.numPositions);
kvyTemp = model.kern.comp{1}.KIC((model.kern.comp{1}.numPositions+1):(2*model.kern.comp{1}.numPositions), ...
    1:model.kern.comp{1}.numPositions);
kyvTemp = model.kern.comp{1}.KIC(1:model.kern.comp{1}.numPositions, ...
    (model.kern.comp{1}.numPositions+1):(2*model.kern.comp{1}.numPositions));
kvvTemp = model.kern.comp{1}.KIC((model.kern.comp{1}.numPositions+1):(2*model.kern.comp{1}.numPositions),  ...
    (model.kern.comp{1}.numPositions+1):(2*model.kern.comp{1}.numPositions));

% We repeat the initial condition matrices, to make easier the process of
% passing them to the sdmultiKernBlock function
howmany = 1 + model.kern.comp{1}.includeVel + model.kern.comp{1}.includeAccel;
kyy = repmat(kyyTemp, howmany, howmany); kvy = repmat(kvyTemp, howmany, howmany);
kyv = repmat(kyvTemp, howmany, howmany); kvv = repmat(kvvTemp, howmany, howmany);

KyuTemp = cell(model.nout,1);

% Compute the Kyy part
for i =1:model.nout,
    covIC(1,1) = kyy(i,i); covIC(2,1) = kvy(i,i);
    covIC(1,2) = kyv(i,i); covIC(2,2) = kvv(i,i);
    switch model.approx
        case {'dtc','fitc', 'dtcvar'}
            model.Kyy{i,1} = sdkernDiagCompute(model.kern.comp{1}.comp{i+1}, model.X{i+1}, covIC);
        case 'pitc'
            model.Kyy{i,1} = real(sdmultiKernComputeBlock(model.kern.comp{1}, model.X{i+1}, i+1, i+1, covIC));
    end
    % Compute the Kyu part
    KyuTemp{i} = multiKernComputeBlock(model.kern.comp{1},  model.X{i+1}, model.X{1}, i+1, 1); 
end

% Compute the Kuu part
KuuTemp = multiKernComputeBlock(model.kern.comp{1}, model.X{1},1,1);

% Organize Kuu and Kyu
for r = 1:model.nlfPerInt,
    model.Kuu{r,1} = KuuTemp{r};
    if isfield(model, 'gamma') && ~isempty(model.gamma)
        model.KuuGamma{r,1} = model.Kuu{r,1} + model.gamma(r)*eye(size(model.Kuu{r,1})); 
    end     
    for i =1:model.nout,
        model.Kyu{i,r} = real(KyuTemp{i}{r});
    end
end


