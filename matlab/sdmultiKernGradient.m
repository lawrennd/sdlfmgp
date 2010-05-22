function g = sdmultiKernGradient(kern, x, x2, covGrad)

% SDMULTIKERNGRADIENT Gradient of SDMULTI kernel's parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the switching
% dynamical multiple output block
% kernel's parameters. As well as the kernel structure and the
% input positions, the user provides a matrix PARTIAL which gives
% the partial derivatives of the function with respect to the
% relevant elements of the kernel matrix.
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG x : the input locations for which the gradients are being
% computed.
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The argument takes
% the form of a square matrix of dimension  numData, where numData is
% the number of rows in X.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters. The ordering of the vector should match
% that provided by the function kernExtractParam.
%
% FORMAT
% DESC computes the derivatives as above, but input locations are
% now provided in two matrices associated with rows and columns of
% the kernel matrix.
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG x1 : the input locations associated with the rows of the
% kernel matrix.
% ARG x2 : the input locations associated with the columns of the
% kernel matrix.
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The matrix should
% have the same number of rows as X1 and the same number of columns
% as X2 has rows.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters.
%
% SEEALSO : multiKernParamInit, kernGradient, multiKernDiagGradient, kernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Pei Gao, 2007
%
% MODIFICATIONS : Mauricio A. Alvarez, 2010

% SDLFMGP

% Divide the matrix of initial conditions in a matrix of proper dimensions
% for positions and velocities

kyyTemp = kern.KIC(1:kern.numPositions, 1:kern.numPositions);
kvyTemp = kern.KIC((kern.numPositions+1):(2*kern.numPositions), 1:kern.numPositions);
kyvTemp = kern.KIC(1:kern.numPositions, (kern.numPositions+1):(2*kern.numPositions));
kvvTemp = kern.KIC((kern.numPositions+1):(2*kern.numPositions),  ...
    (kern.numPositions+1):(2*kern.numPositions));

% We repeat the initial condition matrices, to make easier the process of
% passing them to the sdmultiKernBlock function
howmany = 1 + kern.includeVel + kern.includeAccel;
kyy = repmat(kyyTemp, howmany, howmany); kvy = repmat(kvyTemp, howmany, howmany);
kyv = repmat(kyvTemp, howmany, howmany); kvv = repmat(kvvTemp, howmany, howmany);
covkyy = zeros(howmany*kern.numPositions);
covkyv = zeros(howmany*kern.numPositions);
covkvy = zeros(howmany*kern.numPositions);
covkvv = zeros(howmany*kern.numPositions);


if iscell(x)
    if nargin > 3 && ~iscell(x2)
        error('Time course information is not matched in Cell format!');
    end
    dim1 = zeros(1, kern.numBlocks);
    dim2 = zeros(1, kern.numBlocks);
    arg = cell(1, kern.numBlocks);
    % Collate arguments.
    for i=1:kern.numBlocks
        dim1(i) = size(x{i}, 1);
        arg{i}{1} = x{i};
        if nargin > 3
            dim2(i) = size(x2{i}, 1);
            arg{i}{2} = x2{i};
        else
            dim2(i) = dim1(i);
            arg{i}{2} = arg{i}{1};
            covGrad = x2;
        end
    end
    

    g = zeros(1, size(kern.paramGroups, 1));
    startVal = 1;
    endVal = 0;

    for i = 1:kern.numBlocks
        covIC(1,1) = kyy(i,i); covIC(2,1) = kvy(i,i);
        covIC(1,2) = kyv(i,i); covIC(2,2) = kvv(i,i);
        endVal = endVal + kern.comp{i}.nParams;
        startOne = sum(dim1(1:(i-1)))+1;
        endOne = sum(dim1(1:i));
        startThree = sum(dim2(1:(i-1))) + 1;
        endThree = sum(dim2(1:i));
        if nargin > 3
            [g(1, startVal:endVal), covGradLocal] = sdmultiKernGradientBlock(kern, ...
                arg{i}{:}, covGrad(startOne:endOne, ...
                startThree:endThree), i, i, covIC);
        else
            [g(1, startVal:endVal), covGradLocal] = sdmultiKernGradientBlock(kern, ...
                arg{i}{1}, covGrad(startOne:endOne, ...
                startThree:endThree), i, i, covIC);
        end
        % Assign each local partial derivative to the correct derivative global
        % derivative
        covkyy(i,i) = covGradLocal(1,1);
        covkvy(i,i) = covGradLocal(2,1);
        covkyv(i,i) = covGradLocal(1,2);
        covkvv(i,i) = covGradLocal(2,2);
        startVal2 = 1;
        endVal2 = 0;
        for j = 1:i-1
            covIC(1,1) = kyy(i,j); covIC(2,1) = kvy(i,j);
            covIC(1,2) = kyv(i,j); covIC(2,2) = kvv(i,j);
            endVal2 = endVal2 + kern.comp{j}.nParams;
            if ~isempty(kern.block{i}.cross{j})
                startTwo = sum(dim2(1:(j-1))) + 1;
                endTwo =  sum(dim2(1:j));
                [g1, g2, covGradLocal] = sdmultiKernGradientBlock(kern, arg{i}{1}, ...
                    arg{j}{2}, covGrad(startOne:endOne, ...
                    startTwo:endTwo), i, j, covIC);

                g(1, startVal:endVal) = g(1, startVal:endVal) + 2*g1;
                g(1, startVal2:endVal2) = g(1, startVal2:endVal2) + 2*g2;
            end
            % Assign each local partial derivative to the correct derivative global
            % derivative
            covkyy(i,j) = covGradLocal(1,1);
            covkvy(i,j) = covGradLocal(2,1);
            covkyv(i,j) = covGradLocal(1,2);
            covkvv(i,j) = covGradLocal(2,2);
            covkyy(j,i) = covGradLocal(1,1);
            covkvy(j,i) = covGradLocal(1,2);
            covkyv(j,i) = covGradLocal(2,1);
            covkvv(j,i) = covGradLocal(2,2);
            startVal2 = endVal2 + 1;
        end
        startVal = endVal + 1;
    end

else

    % Collate arguments.
    dim1 = size(x, 1);
    arg{1} = x;
    if nargin > 3
        dim2 = size(x2, 1);
        arg{2} = x2;
    else
        dim2 = dim1;
        covGrad = x2;
    end

    g = zeros(1, size(kern.paramGroups, 1));
    startVal = 1;
    endVal = 0;
    for i = 1:kern.numBlocks
        covIC(1,1) = kyy(i,i); covIC(2,1) = kvy(i,i);
        covIC(1,2) = kyv(i,i); covIC(2,2) = kvv(i,i);        
        endVal = endVal + kern.comp{i}.nParams;
        startOne = (i-1)*dim1 + 1;
        endOne = i*dim1;
        if nargin > 3
            [g(1, startVal:endVal), covGradLocal] = sdmultiKernGradientBlock(kern, ...
                arg{:}, ...
                covGrad(startOne:endOne, ...
                (i-1)*dim2 + 1:i*dim2), ...
                i, i, covIC);
        else
            [g(1, startVal:endVal), covGradLocal] = sdmultiKernGradientBlock(kern, ...
                arg{1}, ...
                covGrad(startOne:endOne, ...
                (i-1)*dim2 + 1:i*dim2), ...
                i, i, covIC);
        end
        % Assign each local partial derivative to the correct derivative global
        % derivative
        covkyy(i,i) = covGradLocal(1,1);
        covkvy(i,i) = covGradLocal(2,1);
        covkyv(i,i) = covGradLocal(1,2);
        covkvv(i,i) = covGradLocal(2,2);
        startVal2 = 1;
        endVal2 = 0;
        for j = 1:i-1
            covIC(1,1) = kyy(i,j); covIC(2,1) = kvy(i,j);
            covIC(1,2) = kyv(i,j); covIC(2,2) = kvv(i,j);
            endVal2 = endVal2 + kern.comp{j}.nParams;
            if ~isempty(kern.block{i}.cross{j})
                startTwo = (j-1)*dim2 + 1;
                endTwo = j*dim2;
                [g1, g2, covGradLocal] = sdmultiKernGradientBlock(kern, ...
                    arg{:}, ...
                    covGrad(startOne:endOne, ...
                    startTwo:endTwo), ...
                    i, j);

                if nargin > 3
                    startThree = (j-1)*dim1 + 1;
                    endThree = j*dim1;
                    [g3, g4, covGradLocal] = sdmultiKernGradientBlock(kern, ...
                        arg{end:-1:1}, ...
                        covGrad(startThree:endThree, ...
                        startTwo:endTwo)', j, i);
                    g(1, startVal:endVal) = g(1, startVal:endVal) + g1 + g4;
                    g(1, startVal2:endVal2) = g(1, startVal2:endVal2) + g2 + g3;
                else
                    g(1, startVal:endVal) = g(1, startVal:endVal) + 2*g1;
                    g(1, startVal2:endVal2) = g(1, startVal2:endVal2) + 2*g2;
                end
            end
            % Assign each local partial derivative to the correct derivative global
            % derivative
            covkyy(i,j) = covGradLocal(1,1);
            covkvy(i,j) = covGradLocal(2,1);
            covkyv(i,j) = covGradLocal(1,2);
            covkvv(i,j) = covGradLocal(2,2);
            covkyy(j,i) = covGradLocal(1,1);
            covkvy(j,i) = covGradLocal(1,2);
            covkyv(j,i) = covGradLocal(2,1);
            covkvv(j,i) = covGradLocal(2,2);
            startVal2 = endVal2 + 1;
        end
        startVal = endVal + 1;
    end

end
% We compute the derivatives with respect to the Cholesky decomposition
finalCovkyy = zeros(kern.numPositions);
finalCovkvy = zeros(kern.numPositions);
finalCovkyv = zeros(kern.numPositions);
finalCovkvv = zeros(kern.numPositions);
startVal1 = 1;
endVal1 = 0;
for i =1:howmany
    endVal1 = endVal1 + kern.numPositions;
    startVal2 = 1;
    endVal2 = 0;
    for j =1:howmany
        endVal2 = endVal2 + kern.numPositions;
        finalCovkyy = finalCovkyy + covkyy(startVal1:endVal1,startVal2:endVal2);
        finalCovkvy = finalCovkvy + covkvy(startVal1:endVal1,startVal2:endVal2);
        finalCovkyv = finalCovkyv + covkyv(startVal1:endVal1,startVal2:endVal2);
        finalCovkvv = finalCovkvv + covkvv(startVal1:endVal1,startVal2:endVal2);
        startVal2 = endVal2 + 1;
    end
    startVal1 = endVal1 + 1;
end
covGradIC = [finalCovkyy finalCovkyv; finalCovkvy finalCovkvv];
Jij = zeros(2*kern.numPositions);
Jji = zeros(2*kern.numPositions);
gIC = zeros(2*kern.numPositions);
for i=1:(2*kern.numPositions)
    for j=1:(2*kern.numPositions)
        Jij(i,j) = 1; Jji(j,i) = 1; 
        gIC(i,j) = sum(sum((kern.LIC*Jij + Jji*kern.LIC').*covGradIC));
        Jij(i,j) = 0; Jji(j,i) = 0;        
    end
end

g = [g(1:kern.nParamsWIC) gIC(:)'];

g = g*kern.paramGroups;
