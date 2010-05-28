function [gKernParam, gX_u] = spsdlfmgpKernGradient(model,dLdKyy, dLdKuy, dLdKuu)

% SPSDLFMGPKERNGRADIENT Gradients of the parameters for the sparse SDLFMGP.
% FORMAT
% DESC computes the gradients of the parameters for sparse switching
% dynamical latent force gp model. When using more than one latent
% function, the partial derivative is organized in a suitable way so that
% the kernGradient routine can be used.
% RETURN gKernParam :  derivatives wrt to the parameters of the kernels of latent,
%    output and independent functions
% RETURN gX_u : derivatives wrt to the inducing points of the latent functions
% ARG model : the model for which the gradients are to be computed.
% ARG dldKyy : the gradient of the likelihood with respect to the block or
% 	   diagonal terms.
% ARG dLdKuu :  the gradient of the likelihood with respect to the
%	   elements of K_uu.
% ARG dLdKuy : the gradient of the likelihood with respect to the
%	   elements of K_uf.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% SDLFMGP

if nargout>1
    gX_u = cell(model.nlf,1);
    for k=1:model.nlf
        gX_u{k} =  zeros(size(model.X{k}));
    end
end

gMulti = zeros(1, size(model.kern.comp{1}.paramGroups, 1));

% Derivatives of the kernel of the latent functions
g1 = real(multiKernGradientBlock(model.kern.comp{1}, model.X{1}, dLdKuu, 1, 1));
gMulti(1:model.kern.comp{1}.comp{1}.nParams) = g1;

% Derivatives of the kernel between latent functions and outputs
startVal = model.kern.comp{1}.comp{1}.nParams + 1;
endVal = model.kern.comp{1}.comp{1}.nParams;
for i = 1:model.nout
    endVal = endVal + model.kern.comp{1}.comp{i+1}.nParams;
    dLdKuyMat = cell2mat(dLdKuy(:,i));
    dLdKyu = mat2cell(dLdKuyMat', model.sizeX(i), model.k);
    g1 = multiKernGradientBlock(model.kern.comp{1}, model.X{i+1}, ...
        model.X{1}, dLdKyu, i+1, 1);
    gMulti(startVal:endVal) = g1;
    startVal = endVal + 1;
end
% % Derivatives of the kernel between ouputs
if ~isempty(dLdKyy{1})
    startVal = model.kern.comp{1}.comp{1}.nParams + 1;
    endVal = model.kern.comp{1}.comp{1}.nParams;
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
    covkyy = zeros(howmany*model.kern.comp{1}.numPositions);
    covkyv = zeros(howmany*model.kern.comp{1}.numPositions);
    covkvy = zeros(howmany*model.kern.comp{1}.numPositions);
    covkvv = zeros(howmany*model.kern.comp{1}.numPositions);
    for i =1:model.nout,
        endVal = endVal + model.kern.comp{1}.comp{i+1}.nParams;
        covIC(1,1) = kyy(i,i); covIC(2,1) = kvy(i,i);
        covIC(1,2) = kyv(i,i); covIC(2,2) = kvv(i,i);
        [g1Out, covGradLocal] = sdmultiKernGradientBlock(model.kern.comp{1}, model.X{i+1}, ...
            dLdKyy{i}, i+1, i+1, covIC);
        gMulti(startVal:endVal) = gMulti(startVal:endVal) + g1Out;
        covkyy(i,i) = covGradLocal(1,1);
        covkvy(i,i) = covGradLocal(2,1);
        covkyv(i,i) = covGradLocal(1,2);
        covkvv(i,i) = covGradLocal(2,2);
        startVal = endVal + 1;
    end
    % We compute the derivatives with respect to the Cholesky decomposition
    finalCovkyy = zeros(model.kern.comp{1}.numPositions);
    finalCovkvy = zeros(model.kern.comp{1}.numPositions);
    finalCovkyv = zeros(model.kern.comp{1}.numPositions);
    finalCovkvv = zeros(model.kern.comp{1}.numPositions);
    startVal1 = 1;
    endVal1 = 0;
    for i =1:howmany
        endVal1 = endVal1 + model.kern.comp{1}.numPositions;
        startVal2 = 1;
        endVal2 = 0;
        for j =1:howmany
            endVal2 = endVal2 + model.kern.comp{1}.numPositions;
            finalCovkyy = finalCovkyy + covkyy(startVal1:endVal1,startVal2:endVal2);
            finalCovkvy = finalCovkvy + covkvy(startVal1:endVal1,startVal2:endVal2);
            finalCovkyv = finalCovkyv + covkyv(startVal1:endVal1,startVal2:endVal2);
            finalCovkvv = finalCovkvv + covkvv(startVal1:endVal1,startVal2:endVal2);
            startVal2 = endVal2 + 1;
        end
        startVal1 = endVal1 + 1;
    end
    covGradIC = [finalCovkyy finalCovkyv; finalCovkvy finalCovkvv];
    Jij = zeros(2*model.kern.comp{1}.numPositions);
    Jji = zeros(2*model.kern.comp{1}.numPositions);
    gIC = zeros(2*model.kern.comp{1}.numPositions);
    for i=1:(2*model.kern.comp{1}.numPositions)
        for j=1:(2*model.kern.comp{1}.numPositions)
            Jij(i,j) = 1; Jji(j,i) = 1;
            gIC(i,j) = sum(sum((model.kern.comp{1}.LIC*Jij + Jji*model.kern.comp{1}.LIC').*covGradIC));
            Jij(i,j) = 0; Jji(j,i) = 0;
        end
    end
else
    gIC = zeros(2*model.kern.comp{1}.numPositions);
end
%gIC = zeros(2*model.kern.comp{1}.numPositions);
gKernParam = [gMulti(1:model.kern.comp{1}.nParamsWIC) gIC(:)'];

if nargout>1
    gX_u = cell2mat(gX_u);
end

