function model = spsdlfmgpUpdateAD(model)

% SPSDLFMGPUPDATEAD Update the representations of A and D associated with 
% the model.
% FORMAT
% DESC updates the representations of A and D in the model when
% called by spmultigpUpdateKernels.
% ARG model : the model for which the representations are being
% updated.
% RETURN model : the model with the A and D representations
% updated.
%
% SEEALSO : spmultigpUpdateKernels, spmultigpExpandParam
%
% COPYRIGHT : Mauricio Alvarez 2010


% SDLFMGP


for r=1:model.nlfPerInt,
    % Kuuinv really is the inverse of Kuu + regularization term
    if isfield(model, 'gamma') && ~isempty(model.gamma)
        [model.Kuuinv{r}, model.sqrtKuu{r}, jitter] = pdinv(model.KuuGamma{r});
        model.logDetKuu{r} = logdet(model.KuuGamma{r}, model.sqrtKuu{r});
    else
        [model.Kuuinv{r}, model.sqrtKuu{r}, jitter] = pdinv(model.Kuu{r});
        model.logDetKuu{r} = logdet(model.Kuu{r}, model.sqrtKuu{r});
    end        
end

for r =1: model.nlfPerInt,
    for k =1: model.nout,
        model.KuuinvKuy{r,k} = model.Kuuinv{r}*model.Kyu{k,r}';
    end
end

for k =1: model.nout,
    KyuKuuinvKuy = zeros(size(model.Kyy{k},1),size(model.Kyy{k},1));
    for r =1: model.nlfPerInt,
        KyuKuuinvKuy = KyuKuuinvKuy + model.Kyu{k,r}*model.KuuinvKuy{r,k};
    end
    switch model.approx
        case {'dtc', 'dtcvar'}
            if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                switch model.noiseOpt
                    case 0
                        model.D{k} = sparseDiag(1/model.beta(k)*ones(size(model.X{k+1},1),1));
                        model.Dinv{k} = sparseDiag(model.beta(k)*ones(size(model.X{k+1},1),1));
                        model.logDetD{k} = -size(model.X{k+1},1)*log(model.beta(k));                  
                    case 1
                        if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
                            model.D{k} = sparseDiag( (1/model.beta(k))*(1./model.nRepeats{k}));
                            model.Dinv{k} = sparseDiag(model.beta(k)*model.nRepeats{k});
                            model.logDetD{k} = sum(log((1/model.beta(k))*(1./model.nRepeats{k})));
                        else
                            error('Model does not contain nRepeats')
                        end                        
                    case 2
                        error('Not implemented yet')
                    case 3
                        % Be aware: the interpretation of beta for this option is as
                        % variance, not precision                       
                        model.D{k} = sparseDiag(1./model.beta{k});
                        model.Dinv{k} = sparseDiag(model.beta{k});
                        model.logDetD{k} = sum(log(1./model.beta{k}));
                end
            else
                model.D{k} = sparseDiag(1/model.beta(k)*ones(size(model.X{k+1},1),1));
                model.Dinv{k} = sparseDiag(model.beta(k)*ones(size(model.X{k+1},1),1));
                model.logDetD{k} = -size(model.X{k+1},1)*log(model.beta(k));
            end
            if strcmp(model.approx, 'dtcvar')
                model.Ktilde{k} = model.Kyy{k} - diag(KyuKuuinvKuy);
            end
        case 'fitc'
            model.D{k} = model.Kyy{k} - diag(KyuKuuinvKuy); % In fitc D is a diagonal matrix.
            if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                switch model.noiseOpt
                    case 0
                        model.D{k} = (model.D{k} + 1/model.beta(k));
                    case 1
                        if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
                            model.D{k} = (model.D{k} + (1/model.beta(k))*(1./model.nRepeats{k}));
                        else
                            error('Model does not contain nRepeats')
                        end
                    case 2
                        error('Not implemented yet')
                    case 3
                        % Be aware: the interpretation of beta for this option is as
                        % variance, not precision
                        model.D{k} = (model.D{k} + 1./model.beta{k});
                end
            else
                model.D{k} = (model.D{k} + 1/model.beta(k));
            end
            model.Dinv{k} = sparseDiag(1./model.D{k});
            model.logDetD{k} = sum(log(model.D{k}));
            model.D{k} = sparseDiag(model.D{k});% This is to keep the flow of the gradients
        case 'pitc'
            model.D{k} = model.Kyy{k} - KyuKuuinvKuy; % In pitc D is a block-diagonal matrix
            if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                switch model.noiseOpt
                    case 0
                         model.D{k} = (model.D{k} + eye(size(model.D{k},1))/model.beta(k));                        
                    case 1
                        if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
                            model.D{k} =(model.D{k} + diag((1/model.beta(k))*(1./model.nRepeats{k})));
                        else
                            error('Model does not contain nRepeats')
                        end
                    case 2
                        error('Not implemented yet')
                    case 3
                        % Be aware: the interpretation of beta for this option is as
                        % variance, not precision
                        model.D{k} = (model.D{k} + diag(1./model.beta{k}));
                end
            else
                model.D{k} = (model.D{k} + eye(size(model.D{k},1))/model.beta(k));
            end
            model.D{k} = checkKernelSymmetry(model.D{k});
            [model.Dinv{k}, model.sqrtD{k}, jitter] = pdinv(model.D{k});
            model.logDetD{k} = logdet(model.D{k}, model.sqrtD{k});
        otherwise
            error('Unknown approximation type')    
    end
    
end

for k =1:model.nlfPerInt,
    for q =1:model.nout,
        switch model.approx
            case {'dtc', 'dtcvar'}
                if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                switch model.noiseOpt
                    case 0
                        model.KuyDinv{k,q} = model.beta(q)*model.Kyu{q,k}';
                    case 1
                        if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
                            model.KuyDinv{k,q} = (diag(model.beta(q)*model.nRepeats{q})*model.Kyu{q,k})';
                        else
                            error('Model does not contain nRepeats')
                        end
                    case 2
                         error('Not implemented yet')
                    case 3           
                         model.KuyDinv{k,q} = (diag(model.beta{q})*model.Kyu{q,k})';
                end
                
                else
                    model.KuyDinv{k,q} = model.beta(q)*model.Kyu{q,k}';
                end                
            case {'fitc','pitc'}
                model.KuyDinv{k,q} = model.Kyu{q,k}'*model.Dinv{q};
        end
    end
end

for r =1:model.nlfPerInt,
    model.KuyDinvy{r,1} = zeros(model.k,1);
    for q =1:model.nout,
        model.KuyDinvy{r} = model.KuyDinvy{r} + model.KuyDinv{r,q}*model.m{q};
    end
end

for k =1:model.nlfPerInt,
    KuyDinvKyu = zeros(model.k);
    for q =1:model.nout,
        KuyDinvKyu = KuyDinvKyu + model.KuyDinv{k,q}*model.Kyu{q,k};
    end        
    if isfield(model, 'gamma') && ~isempty(model.gamma)
        model.A{k,k} = model.KuuGamma{k} + KuyDinvKyu;
    else
        model.A{k,k} = model.Kuu{k} + KuyDinvKyu;
    end
    for r =1:k-1,
        KuyDinvKyu = zeros(model.k, model.k);
        for q =1:model.nout,
            KuyDinvKyu = KuyDinvKyu + model.KuyDinv{k,q}*model.Kyu{q,r};
        end
        model.A{k,r} = KuyDinvKyu;
        model.A{r,k} = KuyDinvKyu';
    end
end

A = cell2mat(model.A);

[Ainv, sqrtA, jitter]  = pdinv(A);
model.logDetA = logdet(A, sqrtA);

model.Ainv  = mat2cell(Ainv, model.k*ones(1,model.nlfPerInt), model.k*ones(1,model.nlfPerInt));
model.sqrtA = mat2cell(sqrtA,model.k*ones(1,model.nlfPerInt), model.k*ones(1,model.nlfPerInt));

%/~
% MAURICIO : Last code
% model.Ainv  = mat2cell(Ainv,model.k*ones(1,model.nlfPerInt), model.k*ones(1,model.nlfPerInt));
% model.sqrtA = mat2cell(sqrtA,model.k*ones(1,model.nlfPerInt),
% model.k*ones(1,model.nlfPerInt)); 
%~/


for r = 1:model.nlfPerInt,
  model.AinvKuyDinvy{r,1} = zeros(model.k,1);
  for k = 1:model.nlfPerInt,
    model.AinvKuyDinvy{r} = model.AinvKuyDinvy{r} + model.Ainv{r,k}*model.KuyDinvy{k};
  end
end


