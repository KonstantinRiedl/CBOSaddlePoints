% CBO-SP numerical example for saddle point problems
%
% This script tests CBO-SP numerically and outputs the approximation to the
% saddle point.
%

%%
clear; clc; close all;

co = set_color();


%% Energy Function E

% % dimension of the ambient space
d1 = 10;
d2 = 10;

% % energy function E
% (E is a function mapping columnwise from R^{d1\times N} \times R^{d2\times N} to R)
objectivefunction = 'SaddleRastrigin';
[E, parametersE, ~, ~] = objective_function(objectivefunction, d1, d2);

% saddle point
xstar = zeros(d1,1);
ystar = zeros(d2,1);


%% Parameters of CBO-SP Algorithm

% time horizon
T = 10;

% discrete time size
dt = 0.01;
 
% number of particles
N = 1000;

% lambda1, lambda2 (parameter of consensus drift term)
lambda1 = 1;
lambda2 = 1;
% type of diffusion
anisotropic = 1;
% sigma1, sigma2 (parameter of exploration term)
sigma1 = sqrt(4);
sigma2 = sqrt(4);

% alpha, beta (weight in Gibbs measure for consensus point computation)
alpha = 10^15;
beta = 10^15;


%% Initialization
X0mean = 1*ones(d1,1);
X0std = 2;
Y0mean = 1*ones(d2,1);
Y0std = 2;


parametersCBOSP = containers.Map({'T', 'dt', 'N', 'alpha', 'beta', 'lambda1', 'lambda2', 'anisotropic', 'sigma1', 'sigma2'},...
                               {  T,   dt,   N,   alpha,   beta,   lambda1,   lambda2,   anisotropic,   sigma1,   sigma2});
parametersInitialization = containers.Map({'X0mean', 'X0std', 'Y0mean', 'Y0std'},...
                                          {  X0mean,   X0std,   Y0mean,   Y0std});


%% CBO Algorithm
%initialization
X0 = X0mean+X0std*randn(d1,N);
X = X0;
Y0 = Y0mean+Y0std*randn(d2,N);
Y = Y0;

% CBO
[xstar_approx, ystar_approx] = CBOSP(E, parametersCBOSP, X0, Y0);

fprintf('saddle point (theoretically): \n')
disp([xstar'])
disp([ystar'])
fprintf('          with objective value: %f\n', E(xstar,ystar))

fprintf('final approximated minimizer  : \n')
disp([xstar_approx'])
disp([ystar_approx'])
fprintf('          with objective value: %f\n', E(xstar_approx,ystar_approx))
if (norm(xstar_approx-xstar) + norm(ystar_approx-ystar))<0.25
    fprintf('************** CBO-SP   successful **************\n')
else
    fprintf('************** CBO-SP UNsuccessful **************\n')
end

