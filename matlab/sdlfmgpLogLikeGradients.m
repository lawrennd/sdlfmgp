function [gParam, gX_u] = sdlfmgpLogLikeGradients(model)

% SDLFMGPLOGLIKEGRADIENTS Compute the gradients for the parameters and X_u.
% FORMAT
% DESC computes the gradients of the sdlfmgp model using as the objective
% function the likelihood (marginal likelihood)
% ARG model : model for which the gradients will be computed
% RETURN gParam : gradients of the parameters
% RETURN gX_u : for sparse approximations, the gradients of the inducing 
% points 
%
% COPYRIGHT : Mauricio A Alvarez, 2010

% SDLFMGP

gX_u = [];
g_scaleBias = gpScaleBiasGradient(model);

switch model.approx
    case 'ftc'
        covGrad = -model.invK + model.invK*model.m*model.m'*model.invK;
        covGrad = 0.5*covGrad;        
        gParam = kernGradient(model.kern, model.X, covGrad);
        if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
             if  strcmp(model.kernType,'lfm') || strcmp(model.kernType,'sim') ...
                    || strcmp(model.kernType,'lfmwhite') ...
                    || strcmp(model.kernType,'simwhite') ...
                    || strcmp(model.kernType,'simnorm') 
                gmuFull = model.m'*model.invK;
                gmu = zeros(1,model.d);
                startVal = 1;
                endVal = 0;
                for j=1:model.d-base,
                    endVal =  endVal + length(model.X{j+base});
                    gmu(j+base) = sum(gmuFull(startVal:endVal));
                    startVal = endVal + 1;
                end
                gmu = gmu(model.nlf+1:end);
                g_meanFunc = meanGradient(model.meanFunction, gmu);
            else
                g_meanFunc = [];
            end
        else
            g_meanFunc = [];
        end
        gParam = [gParam g_scaleBias g_meanFunc];
        if isfield(model, 'fix')
            for i = 1:length(model.fix)
                gParam(model.fix(i).index) = 0;
            end
        end
    case {'dtc', 'fitc', 'pitc', 'dtcvar'}
        % Sparse approximations.
        if isfield(model, 'beta') && ~isempty(model.beta)
            [dLdKyy, dLdKuy, dLdKuu, dLdmu, gBeta] = spmultigpLocalCovGradient(model);
            %[dLdKyy, dLdKuy, dLdKuu, dLdmu, gBeta] = spmultigpLocalCovGradient2(model);
            if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                if model.noiseOpt == 3
                    gBeta = [];
                else
                    fhandle = str2func([model.betaTransform 'Transform']);
                    gBeta = gBeta.*fhandle(model.beta, 'gradfact');
                end
            else
                fhandle = str2func([model.betaTransform 'Transform']);
                gBeta = gBeta.*fhandle(model.beta, 'gradfact');
            end
        else                     
            [dLdKyy, dLdKuy, dLdKuu, dLdmu] = spmultigpLocalCovGradient(model);
            gBeta = [];
        end
        if ~model.fixInducing
            [gParam, gX_u] =  spsdlfmgpKernGradient(model, dLdKyy, dLdKuy, dLdKuu);
        else
            gParam = spsdlfmgpKernGradient(model, dLdKyy, dLdKuy, dLdKuu);
        end
        if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
            if  strcmp(model.kernType,'lfm') || strcmp(model.kernType,'sim') ...
                    || strcmp(model.kernType,'lfmwhite') ...
                    || strcmp(model.kernType,'simwhite') ...
                    || strcmp(model.kernType,'simnorm')
                gmu = zeros(1,model.nout);
                for j=1:model.nout,
                    gmu(j) = sum(dLdmu{j});
                end
                g_meanFunc = meanGradient(model.meanFunction, gmu);
            else
                g_meanFunc = [];
            end
        else
            g_meanFunc = [];
        end
        gParam = [gParam g_scaleBias g_meanFunc gBeta];
        if isfield(model, 'fix')
            for i = 1:length(model.fix)
                gParam(model.fix(i).index) = 0;
            end
        end
        % if there is only one output argument, pack gX_u and gParam into it.
        if ~model.fixInducing || nargout > 1
            gParam = [gParam gX_u(:)'];
        end
end





