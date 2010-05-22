function sdmultiKernDisplay(kern, varargin)

% SDMULTIKERNDISPLAY Display parameters of the SDMULTI kernel.
% FORMAT
% DESC displays the parameters of the switching dynamical multiple output 
% block kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : multiKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS : Mauricio A. Alvarez , 2010
 
% KERN

if nargin > 1
  spacing = repmat(32, 1, varargin{1});
  varargin{1} = varargin{1}+2;
else
  spacing = [];
  varargin{1} = 2;
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Multiple output block kernel:\n')
for i = 1:length(kern.comp)
  fprintf(spacing);
  fprintf('Block %d\n', i)
  kernDisplay(kern.comp{i}, varargin{:});
end
% Print positions covariance
for i = 1:kern.numPositions
    for j = 1:kern.numPositions
        fprintf(spacing)
        fprintf('Covariance at time zero between position %d and position %d: %f\n', ...
            i,j, kern.KIC(i,j));    
    end    
end
% Print velocities against positions covariance
for i = kern.numPositions+1:2*kern.numPositions
    for j = 1:kern.numPositions
        fprintf(spacing)
        fprintf('Covariance at time zero between velocity %d and position %d: %f\n', ...
            i - kern.numPositions, j, kern.KIC(i,j));     
    end    
end
% Print positions against velocities covariance
for i = 1:kern.numPositions
    for j = kern.numPositions+1:2*kern.numPositions
        fprintf(spacing)
        fprintf('Covariance at time zero between position %d and velocity %d: %f\n', ...
            i,j-kern.numPositions, kern.KIC(i,j));     
    end    
end
% Print velocities covariance
for i = kern.numPositions+1:2*kern.numPositions
    for j = kern.numPositions+1:2*kern.numPositions
        fprintf(spacing)
        fprintf('Covariance at time zero between velocity %d and velocity %d: %f\n', ...
            i-kern.numPositions,j-kern.numPositions, kern.KIC(i,j));     
    end    
end


