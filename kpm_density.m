function [F, P] = kpm_density(M, x, y)
%KPM_DENSITY  Evaluate the output of Algorithms 2/3 on a tensor grid.
%
%   [F, P] = kpm_density(M, x, y)
%
% M is the (m+1)-by-(m+1) damped moment matrix returned by sde.m; x and y are
% vectors of points in (-1,1).  The outputs are arranged like meshgrid(x,y),
% i.e. row i corresponds to y(i) and column j to x(j), so they can be passed
% directly to imagesc(x,y,F) or contourf(x,y,F):
%
%   P(i,j) = p_KPM(x_j, y_i) = sum_{a,b} M(a+1,b+1) Ttilde_a(x_j) Ttilde_b(y_i)
%   F(i,j) = P(i,j) / ( sqrt(1-x_j^2) sqrt(1-y_i^2) )
%
% F is the density (with respect to dx dy) of the estimate of mu_A.

    m  = size(M,1) - 1;
    Tx = cheb_evaluation_1d(x, m);           % numel(x)-by-(m+1)
    Ty = cheb_evaluation_1d(y, m);           % numel(y)-by-(m+1)
    P  = Ty * M.' * Tx.';
    F  = P ./ (sqrt(1 - y(:).^2) * sqrt(1 - x(:).^2).');
end
