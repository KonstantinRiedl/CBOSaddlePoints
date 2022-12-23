% Consensus based optimization for saddle point problems (CBO-SP) 
%
% This function performs CBO-SP.
% 
% 
% [xstar_approx, ystar_approx] = CBOSP(E, parametersCBOSP, X0, Y0)
% 
% input:    E                 = objective function E (as anonymous function)
%           parametersCBOSP   = suitable parameters for CBOSP
%                             = [T, dt, N, lambda1, lambda2, anisotropic, sigma1, sigma2, alpha, beta]
%               - T           = time horizon
%               - dt          = time step size
%               - N           = number of particles
%               - lambda1     = consensus drift parameter for minimization
%               - lambda2     = consensus drift parameter for maximization
%               - anisotropic = type of diffusion/noise
%               - sigma1      = exploration/noise parameter for minimization
%               - sigma2      = exploration/noise parameter for maximization
%               - alpha       = weight/temperature parameter alpha for minimization
%               - beta        = weight/temperature parameter beta for maximization
%           X0, Y0            = initial position of the particles
%           
% output:   xstar_approx      = approximation to xstar
%           ystar_approx      = approximation to ystar
%

function [xstar_approx, ystar_approx] = CBOSP(E, parametersCBOSP, X0, Y0)

% get parameters
T = parametersCBOSP('T');
dt = parametersCBOSP('dt');
lambda1 = parametersCBOSP('lambda1');
lambda2 = parametersCBOSP('lambda2');
anisotropic = parametersCBOSP('anisotropic');
sigma1 = parametersCBOSP('sigma1');
sigma2 = parametersCBOSP('sigma2');
alpha = parametersCBOSP('alpha');
beta = parametersCBOSP('beta');

% initialization
X = X0;
Y = Y0;

for k = 1:T/dt
    
    % % CBO-SP iteration (minimization)
    % compute current consensus point x_alpha for minimization
    x_alpha = compute_consensus(E, alpha, X, Y);

    % position updates of one iteration of CBO-SP for minimization
    X = CBOSP_update(dt, lambda1, anisotropic, sigma1, x_alpha, X);

    % % CBO-SP iteration (maximization)
    % compute current consensus point y_alpha for maximization
    y_beta = compute_consensus(E, -beta, X, Y);

    % position updates of one iteration of CBO-SP for maximization
    Y = CBOSP_update(dt, lambda2, anisotropic, sigma2, y_beta, Y);
    
end

xstar_approx = compute_consensus(E, alpha, X, Y);
ystar_approx = compute_consensus(E, -beta, X, Y);

end
