function k = sdkernDiagCompute(kern, x, covIC)

% SDKERNDIAGCOMPUTE Compute the diagonal of switching dynamical LFM kernels
% FORMAT
% DESC computes the diagonal of a kernel matrix for the given kernel.
% ARG kern : kernel structure to be computed.
% ARG X : input data matrix (rows are data points) to the kernel computation.
% RETURN K : vector containing computed diagonal elements of the
% kernel structure.
%
% SEEALSO : kernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2003, 2004, 2005, 2006
%
% MODIFICATIONS : Mauricio A. Alvarez, 2010

% SDLFMGP

fhandle = str2func([kern.type 'KernDiagCompute']);
k = fhandle(kern, x, covIC);
