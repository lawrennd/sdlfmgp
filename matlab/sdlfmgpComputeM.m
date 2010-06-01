function m = sdlfmgpComputeM(model)

% SDLFMGPCOMPUTEM Compute the matrix m of mean values, given the model.
% FORMAT
% DESC computes the matrix m (the scaled, bias and mean function
% removed matrix of the targets), given the model.
% ARG model : the model for which the values are to be computed.
% ARG m : the scaled, bias and mean function removed values.
%
% SEEALSO : sdlfmgpCreate, multigpComputeM
%
% COPYRIGHT : Mauricio Alvarez, 2010

% SDLFMGP

switch model.approx
    case 'ftc'
        if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
            med = meanCompute(model.meanFunction, model.X, model.nlfPerInt);            
            m = model.y - med;
        else
            m = model.y;
        end
        startVal=1;
        endVal=0;
        for j=1:model.d 
            endVal = endVal + size(model.X{j}, 1);
            m(startVal:endVal, 1) = m(startVal:endVal, 1) - model.bias(j);
            if model.scale(j)~=1
                m(startVal:endVal, 1) = m(startVal:endVal, 1)/model.scale(j);
            end
            startVal = endVal+1;
        end
    case {'dtc','fitc','pitc', 'dtcvar'}
        if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)   
            mu = meanCompute(model.meanFunction, mat2cell(ones(model.nout,1), ones(model.nout,1), 1), 0);
        else
            mu = zeros(model.nout,1);
        end
        m = cell(model.nout,1);
        startVal = 1;
        endVal = 0;
        for j=1:model.d
            endVal = endVal + size(model.X{j+1}, 1);
            m{j} = model.y(startVal:endVal) - mu(j);
            m{j} = m{j} - model.bias(j);
            if model.scale(j)~=1
                m{j} = m{j}/model.scale(j);
            end
            startVal = endVal+1;
        end
end


