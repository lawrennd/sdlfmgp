function term  = sdlfmaMeanCompute(sdlfmKern, t , option)

% SDLFMAMEANCOMPUTE Acceleration mean for the switching dynamical LFM model.
% Computes the terms $a_d$ and $b_d$ that appear in the mean function 
% associated with the switching dynamical LFM model. If the mean function 
% is mu(t), then
%
%     mu(t) = a_d(t)y_d(t_0) + b_d(t)\dot{y}_d(t_0),
%
% where $y_d(t_0)$ is the initial condition associated to the position and
% $\dot{y}_d(t_0)$ is the initial condition associated to the velocity.
% 
% FORMAT
% DESC
% ARG sdlfmKern : switching dynamical LFM kernel structure with the
% parameters.
% ARG t : input times for which the mean is to be computed.
% RETURN term : the value of $a_d$.
%
% FORMAT
% DESC
% Computes the terms that appear in the mean function associated with the
% switching dynamical LFM model.
% ARG sdlfmKern : switching dynamical LFM kernel structure with the
% parameters.
% ARG t : input times for which the mean is to be computed.
% ARG option : indicates which term of the mean should be computed. Option
% 'Pos' computes the term $a_d$ and option 'Vel' computes $b_d$ that 
% accompanies the initial condition of the velocity.
% RETURN term : the value of $a_d$ or $b_d$ depending of option.
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
        term = (alpha^2/omega + omega)*exp(-alpha*t).*...
           (alpha*sin(freq) -omega*cos(freq));
    case 'Vel'
        term = exp(-alpha*t).*((alpha^2/omega - omega)*sin(freq)...
            -2*alpha*cos(freq));         
    otherwise
     error('No recognized option')   
end
