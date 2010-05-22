function kyy = organizeIC(kyy, tempPosPos, i, j)

% ORGANIZEIC A helper function that accommodates IC in a specific manner.
% FORMAT
% DESC a helper function that accommodates initial conditions (IC) in a 
% specific manner.
% ARG kyy : current initial condition
% ARG tempPosPos : cell array contaning initial conditions computed for all
% next intervals
% ARG i : current interval in system 1
% ARG j : current interval in system 2
% RETURN kyy : updated initial condition
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% SDLFMGP

kyy(i+1, j+1) = tempPosPos{i,j}(2,2);

if i==1
    kyy(i,j+1) = tempPosPos{i,j}(1,2);
else
    kyy(i+1,j)  = tempPosPos{i,j}(2,1);
end
