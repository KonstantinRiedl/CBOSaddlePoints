% Objective-function function
%
% This function returns the objective function E(x,y) as an anonymous 
% function together with a set of necessary parameters (parametersE) for  
% plotting. Additionally, the function returns two further sets of 
% parameters which are suitable parameters for CBOSP (parametersOptimizer) 
% as well as a suitable initialization (parametersInitialization).
%
% The function E maps columnwise from R^{d1\times N} \times R^{d2\times N} 
% to R^N, i.e., for matrices X in R^{d1\times N} and Y in R^{d2\times N} 
% the function is applied to every column and returns a row vector in R^N 
% (matrix in R^{1\times N}).
% 
% 
% [E, parametersE, parametersOptimizer, parametersInitialization] = objective_function(name, d1, d2)
% 
% input:    name                = name of objective function E
%                                 (e.g., SaddleRastrigin or NonSeparateRastrigin)
%           d1                  = ambient dimension d1 of the minimization
%           d2                  = ambient dimension d2 of the maximization
%           
% output:   E                   = objective function E (as anonymous function)
%           parametersE         = necessary parameters for plotting of E
%                               = [xrange_plot, yrange_plot, zrange_plot]
%           parametersOptimizer = suitable parameters for CBO
%                               = [T, dt, N, lambda, gamma, learning_rate, sigma, alpha]
%               - T             = time horizon
%               - dt            = time step size
%               - N             = number of particles
%               - lambda1       = consensus drift parameter for minimization
%               - lambda2       = consensus drift parameter for maximization
%               - anisotropic   = type of diffusion/noise
%               - sigma1        = exploration/noise parameter for minimization
%               - sigma2        = exploration/noise parameter for maximization
%               - alpha        = weight/temperature parameter alpha for minimization
%               - beta        = weight/temperature parameter alpha for maximization
%           parametersInitialization = suitable parameters for initialization
%                               = [X0mean, X0std, Y0mean, Y0std]
%               - X0mean        = mean of initial distribution for minimzation
%               - X0std         = standard deviation of initial distribution for minimzation
%               - Y0mean        = mean of initial distribution for maximization
%               - Y0std         = standard deviation of initial distribution for maximization
%

function [E, parametersE, parametersOptimizer, parametersInitialization] = objective_function(name, d1, d2)

if nargin == 3
    % standard CBOSP parameters
    if max(d1,d2) == 1
        parametersOptimizer = containers.Map({'T', 'dt', 'N', 'alpha', 'beta', 'lambda1', 'lambda2', 'anisotropic', 'sigma1', 'sigma2'},...
                                             {  4, 0.04, 100,    10^15,    10^15,         1,         1,             1,      0.1,      0.1});
        parametersInitialization = containers.Map({'X0mean', 'X0std', 'Y0mean', 'Y0std'},...
                                                  {       3,       2,        3,       2});
    elseif max(d1,d2) == 2
        parametersOptimizer = containers.Map({'T', 'dt', 'N', 'alpha', 'beta', 'lambda1', 'lambda2', 'anisotropic', 'sigma1', 'sigma2'},...
                                             {  4, 0.02, 100,    10^15,    10^15,         1,         1,             1,      0.1,      0.1});
        parametersInitialization = containers.Map({'X0mean', 'X0std', 'Y0mean', 'Y0std'},...
                                                  {   [4;6],       8,    [4;6],       8});
    else
        parametersOptimizer = [];
        parametersInitialization = [];
    end
else
    error('Input error. Wrong number of inputs.')
end
    
if strcmp(name,'Saddle')
    E = @(x,y) sum(x.*x,1) - sum(y.*y,1);
    xrange_plot = [-2.5;5];
    yrange_plot = [-2.5;5];
    parametersE = [xrange_plot, yrange_plot, [-25;25]];
    parametersOptimizer = parametersOptimizer;
    parametersInitialization = parametersInitialization;
elseif strcmp(name,'SaddleRastrigin')
    E = @(x,y) sum(x.*x,1) + 1.5*sum(1-cos(2*pi*x),1) - (sum(y.*y,1) + 1.5*sum(1-cos(2*pi*y),1));
    xrange_plot = [-2.5;5];
    yrange_plot = [-2.5;5];
    parametersE = [xrange_plot, yrange_plot, [-25;25]];
    parametersOptimizer = parametersOptimizer;
    parametersInitialization = parametersInitialization;
elseif strcmp(name,'SaddleAckley')
    E = @(x,y) -20*exp(-0.2*sqrt(1/d1*sum(x.*x,1)))-exp(1/d1*sum(cos(2*pi*x),1))+exp(1)+20 - (-20*exp(-0.2*sqrt(1/d2*sum(y.*y,1)))-exp(1/d2*sum(cos(2*pi*y),1))+exp(1)+20);
    xrange_plot = [-2.5;5];
    yrange_plot = [-2.5;5];
    parametersE = [xrange_plot, yrange_plot, [-15;15]];
    parametersOptimizer = parametersOptimizer;
    parametersInitialization = parametersInitialization;
elseif strcmp(name,'SaddleRastrigin_x')
    E = @(x,y) sum(x.*x,1) + 1.5*sum(1-cos(2*pi*x),1) - (sum(y.*y,1));
    xrange_plot = [-2.5;5];
    yrange_plot = [-2.5;5];
    parametersE = [xrange_plot, yrange_plot, [-25;25]];
    parametersOptimizer = parametersOptimizer;
    parametersInitialization = parametersInitialization;
    disp('This objective function is not implemented yet.')
elseif strcmp(name,'SaddleRastrigin_Ackley')
    E = @(x,y) sum(x.*x,1) + 1.5*sum(1-cos(2*pi*x),1) - (-20*exp(-0.2*sqrt(1/d2*sum(y.*y,1)))-exp(1/d2*sum(cos(2*pi*y),1))+exp(1)+20);
    xrange_plot = [-2.5;5];
    yrange_plot = [-2.5;5];
    parametersE = [xrange_plot, yrange_plot, [-25;25]];
    parametersOptimizer = parametersOptimizer;
    parametersInitialization = parametersInitialization;
elseif strcmp(name,'NonSeparate')
    E = @(x,y) sum(x.*x,1) - 2*sum(x.*y,1) - sum(y.*y,1);
    xrange_plot = [-2.5;5];
    yrange_plot = [-2.5;5];
    parametersE = [xrange_plot, yrange_plot, [-50;50]];
    parametersOptimizer = parametersOptimizer;
    parametersInitialization = parametersInitialization;
elseif strcmp(name,'NonSeparateRastrigin')
    E = @(x,y) sum(x.*x,1) + 1.5*sum(1-cos(2*pi*x),1) - 2*sum(x.*y,1) - (sum(y.*y,1) + 1.5*sum(1-cos(2*pi*y),1));
    xrange_plot = [-2.5;5];
    yrange_plot = [-2.5;5];
    parametersE = [xrange_plot, yrange_plot, [-50;50]];
    parametersOptimizer = parametersOptimizer;
    parametersInitialization = parametersInitialization;
elseif strcmp(name,'generalSaddle')
    n = 10;
    A = zeros(d1);
    B = zeros(d2,d1);
    C = zeros(d2);
    for i=1:n
        A_i = randn(d1);
        A_i = A_i'*A_i;
        B_i = randn(d2,d1);
        C_i = randn(d2);
        C_i = C_i'*C_i;
        A = A + A_i;
        B = B + B_i;
        C = C + C_i;
    end
    %S = [A, B'; B, -C];
    %E = @(x,y) sum([x;y].*(S*[x;y]),1)/n;
    E = @(x,y) (sum(x.*(A*x),1)/2 + sum(x.*(B'*y),1) - sum(y.*(C*y),1)/2)/n;
    
    xrange_plot = [-10;15];
    yrange_plot = [-10;15];
    parametersE = [xrange_plot, yrange_plot, [-50;50]];
    parametersOptimizer = parametersOptimizer;
    parametersInitialization = parametersInitialization;
else
    disp('This objective function is not implemented yet.')
end