function kern = sdmultiKernExpandParam(kern, params)

% SDMULTIKERNEXPANDPARAM Expands parameters into a SDMULTI kernel struct.
% FORMAT
% DESC returns a multiple output block kernel structure filled with the
% parameters in the given vector. This is used as a helper function to
% enable parameters to be optimised in, for example, the NETLAB
% optimisation functions. It has the same functionality than
% multiKernExpandParam.m, but also expands the parameters associated to the
% initial conditions of the SDLFMGP model kernel.
% ARG kern : the kernel structure in which the parameters are to be
% placed.
% ARG param : vector of parameters which are to be placed in the
% kernel structure.
% RETURN kern : kernel structure with the given parameters in the
% relevant locations.
%
% SEEALSO : multiKernParamInit, multiKernExtractParam, kernExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS : Mauricio A. Alvarez, 2010

% SDLFMGP

params = params*kern.paramGroups';
startVal = 1;
endVal = 0;
for i = 1:length(kern.comp)
  endVal = endVal + kern.comp{i}.nParams;
  if ~isfield(kern, 'fixedBlocks') || ~kern.fixedBlocks(i),
    kern.comp{i} = kernExpandParam(kern.comp{i}, params(1, startVal:endVal));
  end
  startVal = endVal + 1;
end

kern.LIC = reshape(params(startVal:end), 2*kern.numPositions, 2*kern.numPositions);
kern.KIC = kern.LIC*kern.LIC';

