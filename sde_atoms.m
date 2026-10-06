function [Z, w, M] = sde_atoms(A, m, b)
%SDE_ATOMS  Atomic (discretized) spectral density estimate for a normal matrix.
%
%   [Z, w] = sde_atoms(A, m)     uses b ~ Unif(S^{2n-1})  (Algorithm 3)
%   [Z, w] = sde_atoms(A, m, b)  uses the supplied vector b
%
% Runs sde.m with maximum degree m and returns the Gauss-Chebyshev
% discretization of the output measure used in the proof of Theorem 1.1:
% an atomic measure on (m+1)^2 points,
%
%   sum_{i,j} w_i w_j p_KPM(x_i, x_j) delta_{x_i + 1i*x_j},
%
% where x_i, w_i = pi/(m+1) are the nodes and weights of (m+1)-point
% Gauss-Chebyshev quadrature of the first kind.  The quadrature is exact for
% degree <= 2m+1, so the atomic measure has exactly the same Chebyshev moments
% through degree m as the continuous estimate.
%
% Outputs:
%   Z : (m+1)^2-by-1 complex atom locations
%   w : (m+1)^2-by-1 nonnegative masses, sum(w) = 1
%   M : damped moment matrix, as returned by sde.m

    if nargin < 3, b = []; end
    M = sde(A, m, b);

    x  = cos((2*(1:m+1).' - 1) * pi / (2*(m+1)));   % Gauss-Chebyshev nodes
    wq = pi / (m+1);                                % Gauss-Chebyshev weight

    [~, P] = kpm_density(M, x, x);                  % P(i,j) = p_KPM(x_j, x_i)
    [X, Y] = meshgrid(x, x);
    Z = X(:) + 1i*Y(:);
    w = wq^2 * P(:);

    w = max(w, 0);          % p_KPM >= 0 since the Jackson kernel is positive;
    w = w / sum(w);         % explicit normalization in the face of rounding errors
end
