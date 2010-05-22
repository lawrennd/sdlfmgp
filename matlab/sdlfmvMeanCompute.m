function term  = sdlfmvMeanCompute(sdlfmKern, t , option)

% SDLFMVMEANCOMPUTE Velocity mean for the switching dynamical LFM model.
% Computes the terms $g_d$ and $h_d$ that appear in the mean function 
% associated with the switching dynamical LFM model. If the mean function 
% is mu(t), then
%
%     mu(t) = g_d(t)y_d(t_0) + h_d(t)\dot{y}_d(t_0),
%
% where $y_d(t_0)$ is the initial condition associated to the position and
% $\dot{y}_d(t_0)$ is the initial condition associated to the velocity.
% 
% FORMAT
% DESC
% ARG sdlfmKern : switching dynamical LFM kernel structure with the
% parameters.
% ARG t : input times for which the mean is to be computed.
% RETURN term : the value of $g_d$.
%
% FORMAT
% DESC
% Computes the terms that appear in the mean function associated with the
% switching dynamical LFM model.
% ARG sdlfmKern : switching dynamical LFM kernel structure with the
% parameters.
% ARG t : input times for which the mean is to be computed.
% ARG option : indicates which term of the mean should be computed. Option
% 'Pos' computes the term $g_d$ and option 'Vel' computes $h_d$ that 
% accompanies the initial condition of the velocity.
% RETURN term : the value of $g_d$ or $h_d$ depending of option.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% SDLFMGP

if nargin < 3
    option = 'Pos';
end
    
alpha = sdlfmKern.damper/(2*sdlfmKern.mass);
omega = sqrt(sdlfmKern.spring/sdlfmKern.mass-alpha^2);
freq = omega*t;

switch option
    case 'Pos'       
        term = -exp(-alpha*t).*sin(freq)*(alpha^2/omega + omega);
    case 'Vel'
        term = -exp(-alpha*t).*((alpha/omega)*sin(freq)-cos(freq));         
    otherwise
     error('No recognized option')   
end
