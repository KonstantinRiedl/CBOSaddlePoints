% CBO-SP illustrative
%
% This script illustrates the optimization procedure of CBO-SP for 2d
% objective functions.
%

%%
clear; clc; close all;

co = set_color();


%% Settings for Easy Handling and Notes
% 
% decide if time steps require pressing some arbitrary key
manual_steps = 0;
% use pre-set CBO-SP setting (overrides manually chosen parameters)
pre_setparameters = 0;

% 3d plot
spatial_plot = 1;

% save video
savevideo = 0;


%% Energy Function E

% % dimension of the ambient space
d1 = 1;
d2 = 1;

% % energy function E
% (E is a function mapping columnwise from R^{d1\times N}xR^{d2\times N} to R)
objectivefunction = 'SaddleRastrigin';
[E, parametersE, parametersCBOSP, parametersInitialization] = objective_function(objectivefunction, d1, d2);

% range of x (and x and y for plotting)
xrange_plot = parametersE(:,1)';
yrange_plot = parametersE(:,2)';
zrange_plot = parametersE(:,3)';
xrange = 100*xrange_plot;
yrange = 100*yrange_plot;

% saddle point
xstar = zeros(d1,1);
ystar = zeros(d2,1);


%% Parameters of CBO Algorithm

% time horizon
T = 4;

% discrete time size

dt = 0.1;
 
% number of particles
N = 20;

% lambda (parameter of consensus drift term)
lambda1 = 1;
lambda2 = 1;
% type of diffusion
anisotropic = 1;
% sigma (parameter of exploration term)
sigma1 = sqrt(0.1);
sigma2 = sqrt(0.1);

% alpha (weight in Gibbs measure for consensus point computation)
alpha = 10^15;
beta = 10^15;


%% Initialization
X0mean = 2*ones(d1,1);
X0std = 4;
Y0mean = 2*ones(d2,1);
Y0std = 4;


%% Use pre-set setting
if pre_setparameters==1
    T = parametersCBOSP('T');
    dt = parametersCBOSP('dt');
    N = parametersCBOSP('N');
    lambda1 = parametersCBOSP('lambda1');
    lambda2 = parametersCBOSP('lambda2');
    anisotropic = parametersCBOSP('anisotropic');
    sigma1 = parametersCBOSP('sigma1');
    sigma2 = parametersCBOSP('sigma2');
    alpha = parametersCBOSP('alpha');
    beta = parametersCBOSP('beta');
    X0mean = parametersInitialization('X0mean');
    X0std = parametersInitialization('X0std');
    Y0mean = parametersInitialization('Y0mean');
    Y0std = parametersInitialization('Y0std');
else
    parametersCBO = containers.Map({'T', 'dt', 'N', 'alpha', 'beta', 'lambda1', 'lambda2', 'anisotropic', 'sigma1', 'sigma2'},...
                                   {  T,   dt,   N,   alpha,   beta,   lambda1,   lambda2,   anisotropic,   sigma1,   sigma2});
    parametersInitialization = containers.Map({'X0mean', 'X0std', 'Y0mean', 'Y0std'},...
                                              {  X0mean,   X0std,   Y0mean,   Y0std});
end


%% Plotting
scatter_markersize = 75;
plot_markersize = 30;
star_markersize = 13;

% % plot setting
figure('Position', [1200 800 500 400]);
set(gcf,'color','w');
%title('Consensus Based Optimization for Saddlepoint Problems','Interpreter','latex','FontSize',18)  
% % plotting energy function E
[X,Y] = meshgrid(xrange_plot(1):.1:xrange_plot(2),yrange_plot(1):.1:yrange_plot(2));
Z = E(X(:)',Y(:)');
Z = reshape(Z,size(X));

if spatial_plot
    Eplot = surfc(X,Y,Z,'FaceAlpha',0.8);
else
    Eplot = surf(X,Y,Z,'FaceAlpha',0.5);
    Eplot.EdgeColor = 'None';
end
hold on

if ~spatial_plot
    contour(X,Y,Z,100);
end

ax = gca;
ax.FontSize = 14;

if spatial_plot
    view(-25,12.5)
else
    view(2)
end

xlim(xrange_plot)
ylim(yrange_plot)
zlim(zrange_plot)
xticks([-2.5 0 2.5 5])
yticks([-2.5 0 2.5 5])
zticks([-25 0 25])

% way of plotting of all points
if spatial_plot
    F = @(x,y) E(x,y);
else
    F = @(x,y) zrange_plot(2)*ones(size(x)); 
end
%F = @(x) 0*zeros(size(sum(x.*x)));

% % plot global minimizer of energy function E
xystarplot = plot3(xstar, ystar, F(xstar, ystar), '*', 'MarkerSize', star_markersize, 'LineWidth', 1.8, "color", co(5,:));
hold on

%title('Setting','Interpreter','latex','FontSize',18)
if spatial_plot
    legend([xystarplot], 'Saddle point $(x^*,y^*)$','Location','northwest','Interpreter','latex','FontSize',18)
else
    legend([xystarplot], 'Saddle point $(x^*,y^*)$','Location','northwest','Interpreter','latex','FontSize',18)
end

if savevideo
    frame(1) = getframe(gcf);
end
    
pause(dt)
if manual_steps
    pause()
end


%% CBO Algorithm
%initialization
X0 = X0mean+X0std*randn(d1,N);
X = X0;
Y0 = Y0mean+Y0std*randn(d2,N);
Y = Y0;

% plot initial setting
fprintf("t=0\n")
%title(sprintf("CBO-SP at time $t=0$"),'Interpreter','latex','FontSize',18)

XY_plot = scatter3(X0, Y0, F(X0,Y0), scatter_markersize, "MarkerFaceColor", co(3,:), "MarkerEdgeColor", co(3,:));
hold on
if spatial_plot
    legend([xystarplot, XY_plot],'Saddle point $(x^*,y^*)$','Particles $(X_0^i, Y_0^i)$','Location','northwest','Interpreter','latex','FontSize',18)
else
    legend([xystarplot, XY_plot],'Saddle point $(x^*,y^*)$','Particles $(X_0^i, Y_0^i)$','Location','northwest','Interpreter','latex','FontSize',18)
end

if savevideo
    frame(2) = getframe(gcf);
end

%%
% CBO-SP
for k = 1:T/dt
    
    pause(dt)
    if manual_steps
        pause()
        % for exporting snapshots at specific different points in time
        %if k==1 || k==11 || k==21 || k==41 || k==101
        %    pause()
        %end
    end

    t = k*dt;
    %title(sprintf("CBO-SP at time $t=%d$",t),'Interpreter','latex','FontSize',18)
    
    % % CBO-SP iteration (minimization)
    fprintf("Minimization at time t=%d\n", t)
    % compute current consensus point x_alpha
    x_alpha = compute_consensus(E, alpha, X, Y);

    % position updates of one iteration of CBO
    X = CBOSP_update(dt, lambda1, anisotropic, sigma1, x_alpha, X);

    % % Visualization of Minimization Step
    % remove all old plotting objects
    delete(XY_plot)
    if k~=1
        delete(xyalpha_plot)
    end


    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    
    % plotting of particles
    XY_plot = scatter3(X, Y, F(X,Y), scatter_markersize, "MarkerFaceColor", co(3,:), "MarkerEdgeColor", co(3,:));
    % plotting of consensus point
    if k>1
        xyalpha_plot = plot3(x_alpha, y_alpha, F(x_alpha, y_alpha), '.', 'MarkerSize', plot_markersize+5, 'LineWidth', 1.8, "color", co(2,:));

        if spatial_plot
            legend([xystarplot, XY_plot, xyalpha_plot], 'Saddle point $(x^*,y^*)$', 'Particles $(X_t^i, Y_t^i)$', 'Consensus point $(x_{\alpha}(\widehat\rho_{X,t}^N),y_{\beta}(\widehat\rho_{Y,t}^N))$','Location','northwest','Interpreter','latex','FontSize',18)
            %legend([vstarplot, V_plot, valpha_plot], 'Saddle point $(x^*,y^*)$', 'Particles $(X_t^i, Y_t^i)$', 'Consensus point $(x_{\alpha}(\widehat\rho_{X,t}^N),y_{\beta}(\widehat\rho_{Y,t}^N))$','Position',[0.175 0.675 0.1 0.2],'Interpreter','latex','FontSize',18))
        else
            legend([xystarplot, XY_plot, xyalpha_plot], 'Saddle point $(x^*,y^*)$', 'Particles $(X_t^i, Y_t^i)$', 'Consensus point $(x_{\alpha}(\widehat\rho_{X,t}^N),y_{\beta}(\widehat\rho_{Y,t}^N))$','Location','northwest','Interpreter','latex','FontSize',18)
        end
    end
    
    if savevideo
        frame(2*k+1) = getframe(gcf);
    end

    pause(dt)
    if manual_steps
        pause()
    end

    % % CBO-SP iteration (maximization)
    fprintf("Maximization at time t=%d\n", t)
    % compute current consensus point v_alpha
    y_alpha = compute_consensus(E, -beta, X, Y);

    % position updates of one iteration of CBO-SP
    Y = CBOSP_update(dt, lambda2, anisotropic, sigma2, y_alpha, Y);
    
    % Visualization of Maximization Step
    % remove all old plotting objects
    delete(XY_plot)
    if k~=1
        delete(xyalpha_plot)
    end
    
    % plotting of particles
    XY_plot = scatter3(X, Y, F(X,Y), scatter_markersize, "MarkerFaceColor", co(3,:), "MarkerEdgeColor", co(3,:));
    % plotting of consensus point
    xyalpha_plot = plot3(x_alpha, y_alpha, F(x_alpha, y_alpha), '.', 'MarkerSize', plot_markersize+5, 'LineWidth', 1.8, "color", co(2,:));

    if spatial_plot
        legend([xystarplot, XY_plot, xyalpha_plot], 'Saddle point $(x^*,y^*)$', 'Particles $(X_t^i, Y_t^i)$', 'Consensus point $(x_{\alpha}(\widehat\rho_{X,t}^N),y_{\beta}(\widehat\rho_{Y,t}^N))$','Location','northwest','Interpreter','latex','FontSize',18)
        %legend([vstarplot, V_plot, valpha_plot], 'Saddle point $(x^*,y^*)$', 'Particles $(X_t^i, Y_t^i)$', 'Consensus point $(x_{\alpha}(\widehat\rho_{X,t}^N),y_{\beta}(\widehat\rho_{Y,t}^N))$','Position',[0.175 0.675 0.1 0.2],'Interpreter','latex','FontSize',18))
    else
        legend([xystarplot, XY_plot, xyalpha_plot], 'Saddle point $(x^*,y^*)$', 'Particles $(X_t^i, Y_t^i)$', 'Consensus point $(x_{\alpha}(\widehat\rho_{X,t}^N),y_{\beta}(\widehat\rho_{Y,t}^N))$','Location','northwest','Interpreter','latex','FontSize',18)
    end
    
    if savevideo
        frame(2*k+2) = getframe(gcf);
    end
    
end
x_alpha = compute_consensus(E, alpha, X, Y);
y_alpha = compute_consensus(E, -beta, X, Y);
fprintf("global minimizer (numerically): [%d;%d]\n", [xstar, ystar])
fprintf("final consensus point         : [%d;%d]\n", [mean(x_alpha,2), mean(y_alpha,2)])


%% Save Video
if savevideo
    if anisotropic
        video = VideoWriter(['CBOSaddlePoints/images_videos/CBOIllustrative_',objectivefunction,'_anisotropic'],'MPEG-4');
    else
        video = VideoWriter(['CBOSaddlePoints/images_videos/CBOIllustrative_',objectivefunction,'_isotropic'],'MPEG-4');
    end
    video.FrameRate = 8;
    open(video);
    writeVideo(video,frame(1));
    writeVideo(video,frame(2));
    for k = 1:T/dt
        writeVideo(video,frame(k+2));
    end
    close(video);
    % save parameters
    if anisotropic
        save(['CBOSaddlePoints/images_videos/CBOIllustrative_',objectivefunction,'_anisotropic_param'], 'objectivefunction', 'E', 'anisotropic', 'xstar', 'ystar', 'd1', 'd2', 'T', 'dt', 'N', 'alpha', 'beta', 'lambda1', 'lambda2', 'sigma1', 'sigma2', 'X0mean', 'X0std', 'Y0mean', 'Y0std')
    else
        save(['CBOSaddlePoints/images_videos/CBOIllustrative_',objectivefunction,'_isotropic_param'], 'objectivefunction', 'E', 'anisotropic', 'xstar', 'ystar', 'd1', 'd2', 'T', 'dt', 'N', 'alpha', 'beta', 'lambda1', 'lambda2', 'sigma1', 'sigma2', 'X0mean', 'X0std', 'Y0mean', 'Y0std')
    end
end
