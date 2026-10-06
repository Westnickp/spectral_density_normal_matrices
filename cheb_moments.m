function Gamma = cheb_moments(A, b, m)
%CHEB_MOMENTS  Algorithm 1 (ChebMoments): mixed normalized Chebyshev moments.
%
%   Gamma = cheb_moments(A, b, m)
%
% Returns the (m+1)-by-(m+1) matrix of mixed Chebyshev moments of the
% weighted spectral measure mu_{A,b},
%
%   Gamma(j+1,k+1) = int Ttilde_j(x) Ttilde_k(y) dmu_{A,b}(x,y),  0 <= j,k <= m,
%
% where Ttilde_0 = 1/sqrt(pi) and Ttilde_j = sqrt(2/pi)*T_j for j >= 1.
% With H = (A+A')/2 and K = (A-A')/(2i), this is Gamma = X'*Y where
% X = [Ttilde_0(H)b, ..., Ttilde_m(H)b] and Y = [Ttilde_0(K)b, ..., Ttilde_m(K)b].
%
% Inputs:
%   A : n-by-n normal matrix (dense or sparse) with norm(A) <= 1
%   b : n-by-1 starting vector 
%   m : maximum polynomial degree 
%
% Cost: m products with each of H and K, i.e. 4m products with A and A'.
% Storage: two n-by-(m+1) blocks.
%
% A is only touched through products A*v and A'*v, so any matrix with a fast
% matrix-vector product can be substituted by editing applyH/applyK below.

    b  = b(:) / norm(b);
    Ah = A';
    applyH = @(v) (A*v + Ah*v) / 2;
    applyK = @(v) (A*v - Ah*v) / (2i);

    X = cheb_krylov(applyH, b, m);
    Y = cheb_krylov(applyK, b, m);

    % Imaginary part is rounding error. Need to convert from complex to  
    % be real-valued for plotting. 
    Gamma = real(X' * Y);
end

function V = cheb_krylov(applyB, b, m)
%CHEB_KRYLOV  Columns Ttilde_j(B)*b, j = 0..m, by the three-term recurrence.
    V = zeros(numel(b), m+1);
    V(:,1) = b;
    if m >= 1
        V(:,2) = applyB(b);
    end
    for j = 2:m
        V(:,j+1) = 2*applyB(V(:,j)) - V(:,j-1);
    end
    % normalize: Ttilde_0 = 1/sqrt(pi), Ttilde_j = sqrt(2/pi) T_j
    V(:,1)     = V(:,1) / sqrt(pi);
    V(:,2:end) = V(:,2:end) * sqrt(2/pi);
end
