function K = sdmultiKernCompute(kern, varargin)

% SDMULTIKERNCOMPUTE Compute the SDMULTI kernel given the parameters and X.
% FORMAT
% DESC computes the kernel function for the switching dynamical multiple 
% output block kernel given inputs associated with rows and columns. The
% kernel structure must contain the information about initial conditions.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the switching dynamical multiple 
% output block kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : multiKernParamInit, multiKernCompute.m
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

if iscell(varargin{1})
  if length(varargin{1}) ~= kern.numBlocks
    error('Time information is not matched among blocks!');
  end
  dim1 = zeros(1, kern.numBlocks);
  dim2 = zeros(1, kern.numBlocks);
  for i = 1:kern.numBlocks
    dim1(i) = size(varargin{1}{i}, 1);    
    if length(varargin)>1
      if length(varargin{1}) ~= length(varargin{2})
        error('Time information is not matched within the block!');
      end      
      dim2(i) = size(varargin{2}{i}, 1);
    else
      dim2(i) = dim1(i);
    end    
    
  end
    
  K = zeros(sum(dim1), sum(dim2));
  
  for i = 1:kern.numBlocks
    startOne = sum(dim1(1:(i-1))) + 1;
    endOne = sum(dim1(1:i));
    startThree = sum(dim2(1:(i-1))) + 1;
    endThree = sum(dim2(1:i));
    
    covIC(1,1) = kyy(i,i); covIC(2,1) = kvy(i,i);
    covIC(1,2) = kyv(i,i); covIC(2,2) = kvv(i,i);
    
    if length(varargin)<2
      K(startOne:endOne, startThree:endThree) = sdmultiKernComputeBlock(kern, ...
                                                    varargin{1}{i}, i, i, covIC);     
    else
      K(startOne:endOne, startThree:endThree) = sdmultiKernComputeBlock(kern, ...
                                     varargin{1}{i}, varargin{2}{i}, i, i, covIC);      
    end
    for j = 1:i-1
      covIC(1,1) = kyy(i,j); covIC(2,1) = kvy(i,j);
      covIC(1,2) = kyv(i,j); covIC(2,2) = kvv(i,j);
        
      if ~isempty(kern.block{i}.cross{j})
        startTwo = sum(dim2(1:(j-1))) + 1;
        endTwo =  sum(dim2(1:j));
        
        if length(varargin)<2
          K(startOne:endOne, startTwo:endTwo) = sdmultiKernComputeBlock(kern, ...
                                varargin{1}{i}, varargin{1}{j}, i, j, covIC);
          K(startTwo:endTwo, startOne:endOne) = K(startOne:endOne, ...
                                                startTwo:endTwo)';
        else
          K(startOne:endOne, startTwo:endTwo) = sdmultiKernComputeBlock(kern, ...
                                varargin{1}{i}, varargin{2}{j}, i, j, covIC);
          startFour = sum(dim1(1:(j-1))) + 1;
          endFour =  sum(dim1(1:j));
          covIC2(1,1) = kyy(j,i); covIC2(2,1) = kvy(j,i);
          covIC2(1,2) = kyv(j,i); covIC2(2,2) = kvv(j,i);          
          K(startFour:endFour, startThree:endThree) = ...
              sdmultiKernComputeBlock(kern, varargin{2}{i}, varargin{1}{j}, j, i, covIC2)';
        end
      end
    end
  end
 stop = 1;
else
  dim1 = size(varargin{1}, 1);
  
  if length(varargin)>1
    dim2 = size(varargin{2}, 1);
  else
    dim2 = dim1;
  end
  
  K = zeros(kern.numBlocks*dim1, kern.numBlocks*dim2);
  
  for i = 1:kern.numBlocks
    startOne = (i-1)*dim1 + 1;
    endOne = i*dim1;
    startThree = (i-1)*dim2 + 1;
    endThree = i*dim2;
    
    covIC(1,1) = kyy(i,i); covIC(2,1) = kvy(i,i);
    covIC(1,2) = kyv(i,i); covIC(2,2) = kvv(i,i);
    
    K(startOne:endOne, startThree:endThree) = sdmultiKernComputeBlock(kern, ...
                                                      varargin{:}, i, i, covIC);
    for j = 1:i-1
      
      covIC(1,1) = kyy(i,j); covIC(2,1) = kvy(i,j);
      covIC(1,2) = kyv(i,j); covIC(2,2) = kvv(i,j);  
        
      if ~isempty(kern.block{i}.cross{j})
        startTwo = (j-1)*dim2 + 1;
        endTwo = j*dim2;
        K(startOne:endOne, startTwo:endTwo) = sdmultiKernComputeBlock(kern, ...
                                                          varargin{:}, i, j);
        if length(varargin)<2
          K(startTwo:endTwo, startOne:endOne) = K(startOne:endOne, ...
                                                startTwo:endTwo)';
        else
          startFour = (j-1)*dim1 + 1;
          endFour = j*dim1;
          K(startFour:endFour, startThree:endThree) = ...
              sdmultiKernComputeBlock(kern, varargin{end:-1:1}, j, i)';
        end
      end
    end
  end
  
end
