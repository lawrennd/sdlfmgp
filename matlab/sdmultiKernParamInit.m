function kern = sdmultiKernParamInit(kern, options)

% SDMULTIKERNPARAMINIT SDMULTI kernel parameter initialisation.
% The swicthing dynamical multiple output block kernel (SDMULTI) is a 
% wrapper kernel designed to represent the situation where there are 
% several Gaussian processes with correlated outputs and initial conditions
% for the system at time zero. We only change some of the features included
% in the multi Kern structure.
% FORMAT
% DESC initialises the switching dynamical multiple output block kernel 
% structure with some default parameters for the initial conditions
% ARG kern : the kernel structure which requires initialisation.
% ARG options : options for the switching dynamical structure.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, multikernParamInit, sdlfmXsdlfmKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006 
%
% MODIFICATIONS : Mauricio A. Alvarez, 2010

% SDLFMGP

% Reorganize the multi kern structure
kern.comp{1}.type = 'sdmulti';
kern.comp{1}.nParamsWIC = kern.comp{1}.nParams;
kern.comp{1}.includeVel = options.includeVel;
kern.comp{1}.includeAccel = options.includeAccel;
if strcmp(options.approx, 'ftc')
    kern.comp{1}.numPositions = kern.comp{1}.numBlocks/(1 + kern.comp{1}.includeVel + kern.comp{1}.includeAccel);
else
    kern.comp{1}.numPositions = (kern.comp{1}.numBlocks-1)/...
        (1 + kern.comp{1}.includeVel + kern.comp{1}.includeAccel);
end
kern.comp{1}.numVelocities = kern.comp{1}.includeVel*kern.comp{1}.numPositions;
kern.comp{1}.numAccelerations = kern.comp{1}.includeAccel*kern.comp{1}.numPositions;
kern.comp{1}.LIC = eye(2*kern.comp{1}.numPositions);
kern.comp{1}.KIC = kern.comp{1}.LIC*kern.comp{1}.LIC';
kern.comp{1}.nParams = kern.comp{1}.nParams + (2*kern.comp{1}.numPositions)^2;
kern.comp{1}.paramGroups = speye(kern.comp{1}.nParams);

% Reorganize the kern structure: count again the number of parameters
kern.nParams = 0;
for i=1:length(kern.comp)
   kern.nParams = kern.nParams + kern.comp{i}.nParams; 
end

kern.paramGroups = speye(kern.nParams);


