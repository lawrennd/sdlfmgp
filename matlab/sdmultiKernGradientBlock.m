function [g1, g2, covGradLocal] = sdmultiKernGradientBlock(kern, X, X2, covGrad, i, j, covIC)

% SDMULTIKERNGRADIENTBLOCK
% FORMAT
% DESC computes the gradient with respect to parameters for a block
% of a switching dynamical multi-output kernel given two matrices of input. 
% ARG kern : the structure containing the kernel.
% ARG X1 : first set of kernel inputs.
% ARG X2 : second set of kernel inputs.
% ARG covGrad : Gradient of the objective function with respect to
% the relevant portion of the kernel matrix.
% ARG i : the row of the block of the kernel to be computed.
% ARG j : the column of the block of the kernel to be computed.
% ARG covIC : covariance for the initial conditions of the first block
% RETURN g1 : the gradient of the kernel parameters from the first
% kernel in the order provided by the relevant kernExtractParam commands.
% RETURN g2 : the gradient of the kernel parameters from the second
% kernel in the order provided by the relevant kernExtractParam commands.
% RETURN covGradLocal : partial derivatives of the initial conditions for
% the first interval
%
% FORMAT
% DESC compute the gradients of a block of a swicthing dynamical 
% multi-output kernel given a single matrix of input.
% ARG kern : the structure containing the kernel.
% ARG X : first set of kernel inputs.
% ARG covGrad : Gradient of the objective function with respect to
% the relevant portion of the kernel matrix.
% ARG i : the row of the block of the kernel to be computed.
% ARG j : the column of the block of the kernel to be computed.
% ARG covIC : covariance for the initial conditions of the first block
% RETURN g1 : the gradient of the kernel parameters from the first
% kernel in the order provided by the relevant kernExtractParam commands.
% RETURN g2 : the gradient of the kernel parameters from the second
% kernel in the order provided by the relevant kernExtractParam commands.
% RETURN covGradLocal : partial derivatives of the initial conditions for
% the first interval
%
% SEEALSO : multiKernCreate, multiKernGradient, multiKernComputeBlock
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS : Mauricio A. Alvarez, 2010

% SDLFMGP

if nargin < 7
  covIC = j;  
  j = i;
  i = covGrad;
  covGrad = X2;
  X2 = [];
end  

% There is no gradient with respect to fixed blocks.
if isfield(kern, 'fixedBlocks') && kern.fixedBlocks(i) && ...
      kern.fixedBlocks(j),
  g1 = 0;
  g2 = 0;
  return;
end

outArg = 2;
if i == j
  fhandle = [kern.comp{i}.type 'KernGradient'];
  transpose = 0;
  arg{1} = kern.comp{i};
  factors = kernFactors(kern.comp{i}, 'gradfact');
  outArg = 1;
else
  if j<i
    fhandle = [kern.block{i}.cross{j} 'KernGradient'];
    transpose = kern.block{i}.transpose(j);
  else
    fhandle = [kern.block{j}.cross{i} 'KernGradient'];
    transpose = ~kern.block{j}.transpose(i);
  end
  if transpose
    arg{1} = kern.comp{j};
    factors{1} = kernFactors(kern.comp{j}, 'gradfact');
    arg{2} = kern.comp{i};
    factors{2} = kernFactors(kern.comp{i}, 'gradfact');
  else
    arg{1} = kern.comp{i};
    factors{1} = kernFactors(kern.comp{i}, 'gradfact');
    arg{2} = kern.comp{j};
    factors{2} = kernFactors(kern.comp{j}, 'gradfact');
  end
end
fhandle = str2func(fhandle);
arg{end+1} = X;
if ~isempty(X2);
  arg{end+1} = X2;
end
switch outArg
 case 1
  [g1, covGradLocal] = fhandle(arg{:}, covGrad, covIC);
  g1(factors.index) = g1(factors.index).*factors.val;
  g2 = covGradLocal;

 case 2 
  [g1, g2, covGradLocal] = fhandle(arg{:}, covGrad, covIC);
  g1(factors{1}.index) = g1(factors{1}.index).*factors{1}.val;
  g2(factors{2}.index) = g2(factors{2}.index).*factors{2}.val;

  if transpose
    g = g2;
    g2 = g1;
    g1 = g;
  end
 otherwise 
  error('Invalid number of out arguments.')
end