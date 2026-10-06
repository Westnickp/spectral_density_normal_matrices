function T = cheb_evaluation_1d(x, maxdeg)
%CHEB_EVALUATION_1D
% Matrix of normalized Chebyshev polynomials evaluated at points x.
%
% Output T is size (numel(x) x (maxdeg+1)), with
%   T(:,1)   = \tilde T_0(x) = 1/sqrt(pi)
%   T(:,k+1) = \tilde T_k(x) = sqrt(2/pi) * T_k(x), k >= 1

x = x(:);
n = numel(x);

% ordinary Chebyshev values
U = zeros(n, maxdeg+1);
U(:,1) = 1;

if maxdeg >= 1
    U(:,2) = x;
end

for k = 2:maxdeg
    U(:,k+1) = 2*x.*U(:,k) - U(:,k-1);
end

% normalize
T = sqrt(2/pi) * U;
T(:,1) = 1/sqrt(pi);

end