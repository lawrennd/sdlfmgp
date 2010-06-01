function [mu, varsig] = sdlfmgpPosteriorMeanVar(model, X, computeAll)

% SDLFMGPPOSTERIORMEANVAR Mean and variance of the posterior distribution.
% FORMAT
% DESC gives the mean and variance of outputs for the switching dynamical
% latent force model Gaussian process.
% ARG model : the model for which posterior is to be computed.
% ARG X : cell array containing locations where outputs are to be computed.
% RETURN mu : cell array containing mean posterior vectors.
% RETURN varsig : cell array containing the variance posterior vector
%
% DESC gives the mean and variance of outputs for the switching dynamical
% latent force model Gaussian process.
% ARG model : the model for which posterior is to be computed.
% ARG X : cell array containing locations where outputs are to be computed.
% ARG computeAll : a flag to indicate if mean posterior for outputs are to
% be computed. It is true by default.
% RETURN mu : cell array containing mean posterior vectors.
% RETURN varsig : cell array containing the variance posterior vector
%
% COPYRIGHT : Mauricio A Alvarez, 2010
%
% SEEALSO : gpPosteriorMeanVar, multigpPosteriorMeanVar

% SDLFMGP

% If computeAll is true, then compute posterior for latent forces and
% outputs. If computeAll is false, then compute posterior only for latent
% forces. This is introduced particularly for applications in which the
% number of outputs >> number of latent forces and we are only interested
% in latent forces posteriors.

if nargin<3
    computeAll = true;
end

% If X is a vector assume it applies for all outputs.
if ~iscell(X)
    xtemp = X;
    X = cell(1, model.d);
    for i = 1:model.d
        X{i} = xtemp;
    end
end

switch model.approx
    case 'ftc'
        % We need first to create the SDRBF kernel and pass the parameters
        % to it.
        model.sdrbfKern.inverseWidth = model.kern.comp{1}.comp{1}.inverseWidth;
        model.sdrbfKern.switchingTimes = model.kern.comp{1}.comp{1}.switchingTimes;
        mu = cell(model.d+model.nlfPerInt,1);
        varsig = cell(model.d+model.nlfPerInt,1);
        Kfu = zeros(model.N, model.nlfPerInt*length(X{1}));
        startVal = 1;
        endVal = 0;
        for i=1:model.nout
            endVal = endVal + length(model.X{i});
            KTemp = sdlfmXsdrbfKernCompute(model.kern.comp{1}.comp{i}, model.sdrbfKern, model.X{i}, X{i});
            Kfu(startVal:endVal, :) = cell2mat(KTemp);
            startVal = endVal + 1;
        end
        KuuTemp = sdrbfKernCompute(model.sdrbfKern, X{1});
        KuuM = zeros(model.nlfPerInt*length(X{1}));
        startVal = 1;
        endVal = 0;
        for i=1:model.nlfPerInt
            endVal = endVal + length(X{i});
            KuuM(startVal:endVal, startVal:endVal) = KuuTemp{i};
            startVal = endVal + 1;
        end
        muTemp1 = Kfu'*model.alpha;
        KuuPost = KuuM - Kfu'*model.invK*Kfu;
        cont = 0;
        startVal = 1;
        endVal = 0;        
        for i=1:model.nlfPerInt
            cont = cont + 1;
            endVal = endVal + length(X{i});
            mu{cont} = muTemp1(startVal:endVal);
            varsig{cont} = diag(KuuPost);            
        end        
        KX_star = kernCompute(model.kern, model.X, X);
        muTemp = KX_star'*model.alpha;
        diagK = kernDiagCompute(model.kern, X);
        Kinvk = model.invK*KX_star;
        varsigTemp = diagK - sum(KX_star.*Kinvk, 1)';
        if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
            m = meanCompute(model.meanFunction, X, model.nlf);
        end
        startVal=1;
        endVal=0;
        for i=1:length(X)
            endVal = endVal + size(X{i},1);
            if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                mu{i+model.nlfPerInt}  = muTemp(startVal:endVal, 1)*model.scale(i) ...
                    + m(startVal:endVal, 1)+ model.bias(i);
            else
                mu{i+model.nlfPerInt}  = muTemp(startVal:endVal, 1)*model.scale(i) + model.bias(i);
            end
            varsig{i+model.nlfPerInt} = varsigTemp(startVal:endVal, 1)*model.scale(i)*model.scale(i);
            startVal = endVal+1;
        end

    case {'dtc','fitc','pitc','dtcvar'}
        if computeAll
            mu = cell(model.nout+model.nlfPerInt,1);
            varsig = cell(model.nout+model.nlfPerInt,1);
        else
            mu = cell(model.nlfPerInt,1);
            varsig = cell(model.nlfPerInt,1);
        end
        Ku_star_u    = cell(model.nlfPerInt);
        KX_star_X2   = cell(model.nout,model.nlfPerInt);
        KX_star_Diag = cell(model.nout,1);
        % Precomputations
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
            KX_star_Diag{i,1} = sdkernDiagCompute(model.kern.comp{1}.comp{i+1}, X{i}, covIC);
            % Compute the Kyu part
            KyuTemp{i} = multiKernComputeBlock(model.kern.comp{1}, X{i}, model.X{1}, i+1, 1);
        end
        % Compute the Kuu part
        KuuTemp = multiKernComputeBlock(model.kern.comp{1}, X{1}, model.X{1},1,1);
        % Organize Kuu and Kyu
        for r = 1:model.nlfPerInt,
            Ku_star_u{r,1} = KuuTemp{r};
            if isfield(model, 'gamma') && ~isempty(model.gamma)
                Ku_star_u{r,1} =  Ku_star_u{r,1} + model.gamma(r)*eye(size(Ku_star_u{r,1}));
            end
            for i =1:model.nout,
                KX_star_X2{i,r} = real(KyuTemp{i}{r});
            end
        end
        KuuinvAinv = cell(model.nlfPerInt);
        for r =1:model.nlfPerInt,
            for q =1:model.nlfPerInt,
                if r ==q,
                    KuuinvAinv{r,q} = model.Kuuinv{r} - model.Ainv{r,r};
                else
                    KuuinvAinv{r,q} = -model.Ainv{r,q};
                end
            end
        end
        % Posterior for the latent functions
        for j=1:model.nlfPerInt
            mu{j} = Ku_star_u{j}*model.AinvKuyDinvy{j};
            varsig{j} = diag(Ku_star_u{j}*model.Ainv{j,j}*Ku_star_u{j}');% This is because we are only interested in the variances
        end
        % Posterior for the output functions
        if computeAll
            if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                m = meanCompute(model.meanFunction, mat2cell(ones(model.nout,1), ones(model.nout,1), 1), 0);
            end
            DinvKyuAinv = cell(model.nlfPerInt,1);
            for i=1:model.nout,
                DinvKyuAinvKuyDinvy = zeros(size(model.X{i+1},1),1);
                mux_star = zeros(size(X{i},1),1);
                Kx_starXu = zeros(size(X{i},1));
                for r =1:model.nlfPerInt,
                    DinvKyuAinv{r} = zeros(size(model.X{i+1},1),model.k);
                    for q =1:model.nlfPerInt,
                        Kx_starXu = Kx_starXu + KX_star_X2{i,r}*KuuinvAinv{r,q}*KX_star_X2{i,q}';
                        DinvKyuAinv{r} = DinvKyuAinv{r} + model.KuyDinv{q,i}'*model.Ainv{q,r};
                    end
                    mux_star = mux_star + KX_star_X2{i,r}*model.AinvKuyDinvy{r};
                    DinvKyuAinvKuyDinvy = DinvKyuAinvKuyDinvy + model.KuyDinv{r,i}'*model.AinvKuyDinvy{r};
                end
                DinvKyuAinvKuyDinv = zeros(size(model.X{i+1},1));
                DinvKyuAinvKuyStar = zeros(size(model.X{i+1},1), size(X{i},1));
                for r =1:model.nlfPerInt,
                    DinvKyuAinvKuyDinv = DinvKyuAinvKuyDinv +  DinvKyuAinv{r}*model.KuyDinv{r,i};
                    DinvKyuAinvKuyStar = DinvKyuAinvKuyStar +  DinvKyuAinv{r}*KX_star_X2{i,r}';
                end
                if model.includeInd
                    c = length(model.kern.comp);
                    if model.includeNoise
                        c = c - 1;
                    end
                    Kw_star_x2 = real(multiKernComputeBlock(model.kern.comp{c}, X{i}, model.X{i+1}, i+1, i+1));
                    Kw_star_fDinvKyuAinvKuyStar = Kw_star_x2*DinvKyuAinvKuyStar;
                    Kw_star_fKyyinvKfw_star = Kw_star_x2*(model.Dinv{i} - DinvKyuAinvKuyDinv)*Kw_star_x2';
                    covInd = Kw_star_fDinvKyuAinvKuyStar + Kw_star_fDinvKyuAinvKuyStar' + Kw_star_fKyyinvKfw_star;
                    muInd = Kw_star_x2*(model.Dinv{i}*model.m{i} - DinvKyuAinvKuyDinvy);
                    mux_star = mux_star + muInd;
                end
                switch model.kernType
                    case {'gg','ggwhite'}
                        mu{i+model.nlfPerInt} = mux_star*model.scale(i) + model.bias(i);
                    otherwise
                        if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
                            mu{i+model.nlfPerInt} = mux_star*model.scale(i) + m(i) + model.bias(i);
                        else
                            mu{i+model.nlfPerInt} = mux_star*model.scale(i) + model.bias(i);
                        end
                end
                if nargout == 2
                    if model.includeInd
                        if isfield(model, 'beta') && ~isempty(model.beta)
                            varsig{i+model.nlfPerInt} = (KX_star_Diag{i} - diag(Kx_starXu) - diag(covInd) + 1/model.beta(i))*model.scale(i)*model.scale(i);
                        else
                            varsig{i+model.nlfPerInt} = (KX_star_Diag{i} - diag(Kx_starXu) - diag(covInd))*model.scale(i)*model.scale(i);
                        end
                    else
                        if isfield(model, 'beta') && ~isempty(model.beta)
                            if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                                switch model.noiseOpt
                                    case 0
                                        varsig{i+model.nlfPerInt} = (KX_star_Diag{i} - diag(Kx_starXu) + 1/model.beta(i))*model.scale(i)*model.scale(i);
                                    case 1
                                    case 2
                                    case 3
                                        varsig{i+model.nlfPerInt} = (KX_star_Diag{i} - diag(Kx_starXu) )*model.scale(i)*model.scale(i);
                                end                                
                            else
                                varsig{i+model.nlfPerInt} = (KX_star_Diag{i} - diag(Kx_starXu) + 1/model.beta(i))*model.scale(i)*model.scale(i);
                            end
                        else
                            varsig{i+model.nlfPerInt} = (KX_star_Diag{i} - diag(Kx_starXu))*model.scale(i)*model.scale(i);
                        end
                    end
                end
            end
        end
    otherwise
        error('multigpPosteriorMeanVar not yet implemented');
end




