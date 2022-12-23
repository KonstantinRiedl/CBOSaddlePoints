% Objective function plot
%
% This script produces plots of the objective function E together with its
% global saddlepoint.
%

%%
clear; clc; close all;

co = set_color();


%% Settings for Easy Handling and Notes
% save plot
pdfexport = 0;


%% Energy Function 

% % dimension of the ambient space
d1 = 1;
d2 = 1;

% % energy function E
% (E is a function mapping columnwise from R^{d1\times N} \times R^{d2\times N} to R)
objectivefunction = 'NonSeparateRastrigin';
[E, parametersE, parametersCBO, parametersInitialization] = objective_function(objectivefunction, d1, d2);

% range of x (and x and y for plotting)
xrange_plot = parametersE(:,1)';
yrange_plot = parametersE(:,2)';
zrange_plot = parametersE(:,3)';
xrange = 100*xrange_plot;
yrange = 100*yrange_plot;

% saddle point
xstar = zeros(d1,1);
ystar = zeros(d2,1);

%% Plotting

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% % plot setting
f = figure('Position', [1200 800 500 400]);

% % plotting energy function E
[X,Y] = meshgrid(xrange_plot(1):.1:xrange_plot(2),yrange_plot(1):.1:yrange_plot(2));
Z = E(X(:)',Y(:)');
Z = reshape(Z,size(X));

%Eplot = surf(X,Y,Z,'FaceAlpha',0.55); % 0.5 und 0.25
%Eplot.EdgeColor = 'None';
%hold on
%contour(X,Y,Z,20);
Eplot = surfc(X,Y,Z,'FaceAlpha',0.55,'EdgeColor','interp'); % 0.5 und 0.25
hold on

view(-25,12.5)
xlim(xrange_plot)
ylim(yrange_plot)
zlim(zrange_plot)

xticks([-2.5 0 2.5 5])
yticks([-2.5 0 2.5 5])
zticks([-50 -25 0 25 50])


% % plot saddle point of energy function E
xystarplot = plot3(xstar, ystar, E(xstar, ystar), '*', 'MarkerSize', 10, 'LineWidth', 1.8, "color", co(5,:));
hold on

ax = gca;
ax.FontSize = 11;


%% Save Image
if pdfexport
    
    %cleanfigure;
    %matlab2tikz('myfile.tex');
    
    print(f,['CBOSaddlePoints/images_videos/ObjectiveFunction2d_',objectivefunction],'-dpdf');

    % save parameters
    save(['CBOSaddlePoints/images_videos/ObjectiveFunction2d_',objectivefunction,'_param'], 'objectivefunction', 'E', 'xstar', 'ystar', 'd1', 'd2')

end

