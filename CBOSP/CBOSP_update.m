% Position updates of one iteration of consensus based optimization for
% saddle point problems (CBO-SP)
%
% This function performs the position updates of one iteration of CBO-SP.
% 
% 
% [XY] = CBOSP_update(E, lambda, anisotropic, sigma, xy_alpha, XY)
% 
% input:    E             = objective function E (as anonymous function)
%           lambda        = consensus drift parameter
%           anisotropic   = type of noise
%           sigma         = exploration/noise parameter
%           xy_alpha      = current empirical consensus point
%           XY            = former x or y positions of the particles
%           
% output:   XY            = x or y positions of the particles afterwards
%

function [XY] = CBOSP_update(dt, lambda, anisotropic, sigma, xy_alpha, XY)

% get parameters
d = size(XY,1);
N = size(XY,2);


% Brownian motion for exploration term
dB = randn(d,N);


% % particle update step (according to SDE)
% consensus drift term
XY = XY - lambda*(XY-xy_alpha)*dt;
% exploration/noise term
if anisotropic
    XY = XY + sigma*abs(XY-xy_alpha)*sqrt(dt).*dB;
else
    XY = XY + sigma*vecnorm(XY-xy_alpha,2,1)*sqrt(dt).*dB;
end

end
