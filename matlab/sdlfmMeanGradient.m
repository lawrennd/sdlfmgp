function [gradAlpha, gradOmega] = sdlfmMeanGradient(sdlfmKern, t , option)

% SDLFMMEANGRADIENT Gradients wrt parameters of the Position mean SDLFM.
% FORMAT
% DESC computes the gradients of the terms $c_d$ that appear in the 
% mean function with respect to $\alpha_d$ and $\omega_d$. 
% ARG sdlfmKern : switching dynamical LFM kernel structure with the
% parameters.
% ARG t : input times for which the mean is to be computed.
% RETURN gradAlpha : gradient of $c_d$ wrt $\alpha_d$
% RETURN gradOmega : gradient of $c_d$ wrt $\omega_d$
%
% FORMAT
% DESC computes the gradients of the terms $c_d$ and $e_d$ that appear in 
% the mean function with respect to $\alpha_d$ and $\omega_d$. 
% ARG sdlfmKern : switching dynamical LFM kernel structure with the
% parameters.
% ARG t : input times for which the mean is to be computed.
% ARG option : indicates which to which term of the mean should be should
% the derivatives wrt $\alpha_d$ and  $\omega_d$ should be computed. Option
% 'Pos' computes the respective derivatives of the term $c_d$ and option 
% 'Vel' computes the respective derivatives of $e_d$.
% RETURN gradAlpha : gradient of $c_d$ or $e_d$ wrt $\alpha_d$
% RETURN gradOmega : gradient of $c_d$ or $e_d$ wrt $\omega_d$
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% SDLFMGP

if nargin < 3
    option = 'Pos';
end
    
alpha = sdlfmKern.damper/(2*sdlfmKern.mass);
omega = sqrt(sdlfmKern.spring/sdlfmKern.mass-alpha^2);
freq = omega*t;
ed = sdlfmMeanCompute(sdlfmKern, t, 'Vel');

switch option
    case 'Pos'
        cd = sdlfmMeanCompute(sdlfmKern, t);
        gradAlpha = -t.*cd + ed;
        gradOmega = exp(-alpha*t).*(t.*((alpha/omega)*cos(freq)-sin(freq)) ...
            - (alpha/omega^2)*sin(freq));      
    case 'Vel'
        gradAlpha = -t.*ed;
        gradOmega = (1/omega)*(t.*exp(-alpha*t).*cos(freq) - ed);
    otherwise
     error('No recognized option')   
end
