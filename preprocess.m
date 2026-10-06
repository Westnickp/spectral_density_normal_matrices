function [c1, c2, nrmb_est] = preprocess(A, k, eta)
%PREPROCESS  Estimate shift/scale so c1*(A-c2*I) has norm < 1.
%
%   [c1, c2, nrmb_est] = preprocess(A, k, eta)
%
% Inputs:
%   A   : square matrix
%   k   : number of power iterations on (A-c2*I)'*(A-c2*I)
%   eta : safety factor in (0,1), default 0.95
%
% Outputs:
%   c2       : estimated center (trace(A)/n via one random probe)
%   nrmb_est : estimated ||A - c2 I||_2
%   c1       : scaling = eta / nrmb_est

if nargin < 3 || isempty(eta)
    eta = 0.95;
end

n = size(A,1);

% random probe for trace estimate
b = randn(n,1) + 1i*randn(n,1);
b = b / norm(b);

% estimate trace(A)/n
c2 = b' * (A*b);

% centered matrix
B = A - c2*speye(n);

% power iteration on B'*B
v = randn(n,1) + 1i*randn(n,1);
v = v / norm(v);

for j = 1:k
    v = B' * (B * v);
    nv = norm(v);
    if nv == 0
        nrmb_est = 0;
        c1 = 1;
        return
    end
    v = v / nv;
end

nrmb_est = norm(B*v);

if nrmb_est == 0
    c1 = 1;
else
    c1 = eta / nrmb_est;
end

end