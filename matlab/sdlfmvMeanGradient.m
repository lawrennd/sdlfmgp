function [gradAlpha, gradOmega] = sdlfmvMeanGradient(sdlfmKern, t , option)

% SDLFMVMEANGRADIENT Gradients wrt parameters of the velocity mean SDLFM.
% FORMAT
% DESC computes the gradients of the terms $g_d$ that appear in the 
% mean function with respect to $\alpha_d$ and $\omega_d$. 
% ARG sdlfmKern : switching dynamical LFM kernel structure with the
% parameters.
% ARG t : input times for which the mean is to be computed.
% RETURN gradAlpha : gradient of $g_d$ wrt $\alpha_d$
% RETURN gradOmega : gradient of $g_d$ wrt $\omega_d$
%
% FORMAT
% DESC computes the gradients of the terms $g_d$ and $h_d$ that appear in 
% the mean function with respect to $\alpha_d$ and $\omega_d$. 
% ARG sdlfmKern : switching dynamical LFM kernel structure with the
% parameters.
% ARG t : input times for which the mean is to be computed.
% ARG option : indicates which to which term of the mean should be should
% the derivatives wrt $\alpha_d$ and  $\omega_d$ should be computed. Option
% 'Pos' computes the respective derivatives of the term $g_d$ and option 
% 'Vel' computes the respective derivatives of $h_d$.
% RETURN gradAlpha : gradient of $g_d$ or $h_d$ wrt $\alpha_d$
% RETURN gradOmega : gradient of $g_d$ or $h_d$ wrt $\omega_d$
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
        gd = sdlfmvMeanCompute(sdlfmKern, t);
        gradAlpha = -t.*gd - 2*alpha*ed;
        gradOmega = -t.*exp(-alpha*t).*cos(freq)*(alpha^2/omega + omega) ...
            - exp(-alpha*t).*sin(freq)*(1- alpha^2/omega^2);
    case 'Vel'
        hd = sdlfmvMeanCompute(sdlfmKern, t, 'Vel');
        gradAlpha = -t.*hd - ed;
        gradOmega = -exp(-alpha*t).*(t.*((alpha/omega)*cos(freq)+sin(freq)) ...
            - (alpha/omega^2)*sin(freq));
    otherwise
     error('No recognized option')   
end
