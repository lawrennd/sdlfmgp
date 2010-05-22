function options = sdlfmgpOptions(approx)

% SDLFMGPOPTIONS Return default options for the SDLFMGP examples.
% FORMAT
% DESC returns the default options in a structure for a SDLFMGP model.
% ARG approx : approximation type, either 'none' (no approximation),
% 'fitc' (fully independent training conditional), 'pitc' (partially
% independent training conditional or 'dtcvar' (variational deterministic
% training conditional).
% RETURN options : structure containing the default options for the
% given approximation type.
%
% SEEALSO : sdlfmgpCreate
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% SDLFMGP

  options.type = 'sdlfmgp';
  
  if nargin<1
    options.approx = 'ftc';
  else
    options.approx = approx;
  end
  % Basic kernel
  options.kernType = 'sdlfm';
  % No include velocity
  options.includeVel = false;
  % No include acceleration
  options.includeAccel = false;
  % Include noise in the model
  options.includeNoise = true;
  %Learn the scales
  options.learnScales = false;
  % Include independent kernel in the model
  options.includeInd = false;
  % Include options to tie the parameters
  options.tieOptions = multigpTieOptions;
  % Method for optimization
  options.optimiser = 'scg';
  % One latent function per interval
  options.nlfPerInt = 1;
  % Set the number of intervals
  options.nIntervals = 3;
  
  % Set options for the kernel
  options.kern.nIntervals = options.nIntervals;
  options.kern.nlfPerInt = options.nlfPerInt;
  options.kern.switchingTimes = [-0.1 0.5 1];
  options.kern.includeVel = false;
  options.kern.includeAccel = false;
  
  % Set to a given mean function to have a mean function.
  options.meanFunction = [];
  % Options structure for mean function options.
  options.meanFunctionOptions = [];
  
  
  
  switch options.approx
   case 'ftc'
    options.numActive = [];
   case {'dtc','fitc','pitc', 'dtcvar'}
    options.numActive = 15;
    options.fixInducing = true;
    options.fixIndices = 1:options.numActive;
    options.includeScalesp = 0;
    if (strcmp(options.approx, 'dtc') || strcmp(options.approx, 'dtcvar'))
        options.beta = 1e-3;
    else
        options.beta = 1e3;
    end
    % Initial position of the inducing variables. Options are 'random',
    % in which random initial locations taken from the used data;
    % 'espaced' the initial locations are equally spaced chosen across
    % all dimensions; 'fixIndices' the indices in the fixIndices
    % options are employed; 'kmeans' the kmeans method gives the
    % initial positions.        
    options.initialInducingPositionMethod = 'espacedInRange'; 
    if strcmp(options.approx, 'dtcvar')
        % Learns the sensitivities variationally
        options.varS = false;
    end        
  end
end


function options = multigpTieOptions
  options.tieIndices = false;
  options.selectMethod = 'free';
end

