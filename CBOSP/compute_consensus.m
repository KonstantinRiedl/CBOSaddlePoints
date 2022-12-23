% Consensus point computation
%
% This function computes the current consensus point based on the positions
% of a given number N of particles xy_i according to the formula
%
%    xy_alpha = sum_i=1^N xy_i*w_alpha(xy_i)/(sum_j=1^N w_alpha(xy_j)),
%
% where w_alpha(xy_i) = exp(-alpha*E_avg(x_i)) with E_avg being
%     E_avg(x) = E(x,y_average) in case of minimization (i.e., alpha>0) or
%     E_avg(y) = E(x_average,y) in case of maximization (i.e., alpha<0).
% Notice that for alpha>0 the function computes the consensus point for a
% minimization problem while for alpha<0 it is computed for a maximization
% problem.
% The implementation is moreover (numerically) stabilized.
% 
% 
% [xy_alpha] = compute_consensus(E, alpha, XY, YX)
% 
% input:    E             = objective function E (as anonymous function)
%           alphabeta     = weight/temperature parameter alpha
%           X             = positions x_i of particles used for computation
%                           of current empirical consensus point xy_alpha
%           Y             = positions y_i of particles used for computation
%                           of current empirical consensus point xy_alpha
%           
% output:   xy_alphabeta  = current empirical consensus point
%

function [xy_alphabeta] = compute_consensus(E, alphabeta, X, Y)

% computation of current empirical consensus point xy_alpha
if alphabeta>0
    Es = E(X, mean(Y,2));
    Emin = min(Es);
    w_alphabeta = exp(-alphabeta*(Es-Emin));
    xy_alphabeta = sum((X.*w_alphabeta),2);
else
    Es = E(mean(X,2), Y);
    Emax = max(Es);
    w_alphabeta = exp(-alphabeta*(Es-Emax));
    xy_alphabeta = sum((Y.*w_alphabeta),2);
end
xy_alphabeta = 1/sum(w_alphabeta)*xy_alphabeta;

end
