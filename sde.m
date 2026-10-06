function [M, Gamma, b] = sde(A, m, b)
%SDE  Algorithm 3 (SpecDens): spectral density estimate for a normal matrix.
%
%   [M, Gamma, b] = sde(A, m)      uses b ~ Unif(S^{2n-1})  (Algorithm 3)
%   [M, Gamma, b] = sde(A, m, b)   uses the supplied vector b (Algorithm 2)
%
% Inputs:
%   A : n-by-n normal matrix with norm(A) <= 1 
%   m : maximum polynomial degree (as in the paper)
%   b : optional starting vector (normalized internally)
%
% Outputs:
%   M     : (m+1)-by-(m+1) damped moments, M(j+1,k+1) = rho_j rho_k Gamma_jk.
%           These are the coefficients of the KPM polynomial
%               p_KPM(x,y) = sum_{j,k} M(j+1,k+1) Ttilde_j(x) Ttilde_k(y),
%           and the estimate of mu_A is the measure with density
%               p_KPM(x,y) / ( sqrt(1-x^2) sqrt(1-y^2) )   on [-1,1]^2,
%           identifying x + iy with a point of the complex plane.
%           Evaluate it with kpm_density.m, or discretize it with sde_atoms.m.
%   Gamma : undamped mixed Chebyshev moments of mu_{A,b} (Algorithm 1)
%   b     : the (normalized) starting vector that was used
%
% Cost: 4m matrix-vector products with A and A'.

    n = size(A,1);
    if nargin < 3 || isempty(b)
        b = randn(n,1) + 1i*randn(n,1);      % uniform direction on S^{2n-1}
    end
    b = b(:) / norm(b);

    Gamma = cheb_moments(A, b, m);
    rho   = jackson_coeffs(m);
    M     = Gamma .* (rho * rho.');
end
