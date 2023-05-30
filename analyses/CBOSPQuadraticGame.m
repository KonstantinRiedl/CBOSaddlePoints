%% Energy Function E

% % dimension of the ambient space
d1 = 40;
d2 = 20;

% % energy function E
% (E is a function mapping columnwise from R^{d1\times N} \times R^{d2\times N} to R)
objectivefunction = 'generalSaddle';
[E, parametersE, ~, ~] = objective_function(objectivefunction, d1, d2);

% saddle point
xstar = zeros(d1,1);
ystar = zeros(d2,1);


%% Parameters of CBO-SP Algorithm

% time horizon
T = 100;

% discrete time size
dt = 0.1;
 
% number of particles
N = 120;

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
X0mean = 4*ones(d1,1);
X0std = 2;
Y0mean = 4*ones(d2,1);
Y0std = 2;


parametersCBOSP = containers.Map({'T', 'dt', 'N', 'alpha', 'beta', 'lambda1', 'lambda2', 'anisotropic', 'sigma1', 'sigma2'},...
                                 {  T,   dt,   N,   alpha,   beta,   lambda1,   lambda2,   anisotropic,   sigma1,   sigma2});
parametersInitialization = containers.Map({'X0mean', 'X0std', 'Y0mean', 'Y0std'},...
                                          {  X0mean,   X0std,   Y0mean,   Y0std});


%% CBO Algorithm
%initialization
M = 10;
success_count = 0;
avg_error = 0;
avg_runtime = 0;

for m = 1:M
    X0 = X0mean+X0std*randn(d1,N);
    X = X0;
    Y0 = Y0mean+Y0std*randn(d2,N);
    Y = Y0;
    
    % CBO
    tic;
    [xstar_approx, ystar_approx] = CBOSP(E, parametersCBOSP, X0, Y0);
    runtime = toc;
    
    fprintf('saddle point (theoretically): \n')
    disp([xstar'])
    disp([ystar'])
    fprintf('          with objective value: %f\n', E(xstar,ystar))
    
    fprintf('final approximated minimizer  : \n')
    disp([xstar_approx'])
    disp([ystar_approx'])
    fprintf('          with objective value: %f\n', E(xstar_approx,ystar_approx))
    %if (norm(xstar_approx-xstar) + norm(ystar_approx-ystar))<0.25
    if (max(norm(xstar_approx-xstar,"inf"), norm(ystar_approx-ystar,"inf")))<10^-3
        fprintf('************** CBO-SP   successful **************\n')
    else
        fprintf('************** CBO-SP UNsuccessful **************\n')
    end
    
    success_count = success_count + ((max(norm(xstar_approx-xstar,"inf"), norm(ystar_approx-ystar,"inf")))<10^-3);
    avg_error = avg_error + max(norm(xstar_approx-xstar,"inf"), norm(ystar_approx-ystar,"inf"));
    avg_runtime = avg_runtime + runtime;

end

1000*avg_runtime/M
avg_error/M
100*success_count/M

