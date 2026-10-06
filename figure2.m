%FIGURE2  Reproduces Figure 2: earth mover's (Wasserstein-1) error of
% Algorithm 3 as a function of the maximum degree m = 2, 4, ..., 128.
%
% W1 is computed exactly by w1_pot.m, which
% calls the network simplex solver in the Python package POT:
%
%   R. Flamary et al., POT: Python Optimal Transport, JMLR 2021
%
% TEST MATRIX.  A = diag(z), n = 1000.  55% of the eigenvalues are equispaced
% on the circle of radius 0.95 centred at 0; the remaining 45% are split
% equally among three Gaussian clusters with centres -0.5, -0.5i, 0.6+0.3i,
% whose real and imaginary parts are independent with standard deviation 0.1
% (conditioned on |z| <= 0.95, so that norm(A) <= 1).
%
% CURVES (all on the same matrix):
%  exact b   b = ones(n,1)/sqrt(n), so mu_{A,b} = mu_A exactly and only the
%            O(1/m) kernel-polynomial error remains.
%  random b  Algorithm 3 as stated, with b ~ Unif(S^{2n-1}).  Tracks the exact
%            curve until it reaches the floor W1(mu_A, mu_{A,b}), then plateaus.
%  floor     W1(mu_A, mu_{A,b}) for the same random b, which does not
%            depend on m.
%
% The scored measure is the (m+1)^2-atom Gauss-Chebyshev discretization of the
% estimate (sde_atoms.m), which has the same Chebyshev moments through degree
% m.
%
% Requires: sde_atoms.m, sde.m, cheb_moments.m, jackson_coeffs.m,
%           cheb_evaluation_1d.m, kpm_density.m, w1_pot.m
% NOTE:     Additionally requires Python with POT installed (pip install pot), configured in
%           MATLAB with pyenv. See README.txt for details on this. 
% Output:   figure2.png, figure2_data.mat

clear; close all

MS = 2.^(1:7);                 % m = 2, 4, ..., 128
N  = 1000;                     % matrix dimension

%% ---------------------------------------------------------- test matrix ----
rng(1);
nb = round(0.15*N);  nc = N - 3*nb;
centres = [-0.5; -0.5i; 0.6+0.3i];
z = zeros(0,1);
for j = 1:3
    % Gaussian cluster, conditioned on lying in the disk |z| <= 0.95 so that
    % norm(A) <= 1 (rejection sampling; affects only a handful of draws)
    w = zeros(0,1);
    while numel(w) < nb
        c = centres(j) + 0.1*(randn(nb,1) + 1i*randn(nb,1));
        w = [w; c(abs(c) <= 0.95)]; 
    end
    z = [z; w(1:nb)]; 
end
z = [z; 0.95*exp(2i*pi*(0:nc-1).'/nc)];     % equispaced circle points
A = spdiags(z, 0, N, N);
a = ones(N,1) / N;             % masses of mu_A

tall = tic;

%% ------------------------------------------------------- (a) exact b ------
bex = ones(N,1) / sqrt(N);
Wex = zeros(size(MS));  Tex = zeros(size(MS));
for i = 1:numel(MS)
    t = tic;
    [Z, w] = sde_atoms(A, MS(i), bex);
    Wex(i) = w1_pot(z, a, Z, w);
    Tex(i) = toc(t);
    fprintf('exact   m = %4d  atoms = %6d  W1 = %.7f  m*W1 = %.3f  (%.1fs)\n', ...
            MS(i), numel(w), Wex(i), MS(i)*Wex(i), Tex(i));
end

%% ------------------------------------------------------ (b) random b ------
rng(2);
b = randn(N,1) + 1i*randn(N,1);  b = b / norm(b);

% the floor W1(mu_A, mu_{A,b}): both measures sit on the same N atoms
t = tic;
Wfloor = w1_pot(z, a, z, abs(b).^2);
fprintf('random  floor W1(mu_A, mu_{A,b}) = %.7f  (%.1fs)\n', Wfloor, toc(t));

Wrand = zeros(size(MS));  Trand = zeros(size(MS));
for i = 1:numel(MS)
    t = tic;
    [Z, w] = sde_atoms(A, MS(i), b);
    Wrand(i) = w1_pot(z, a, Z, w);
    Trand(i) = toc(t);
    fprintf('random  m = %4d  atoms = %6d  W1 = %.7f  (%.1fs)\n', ...
            MS(i), numel(w), Wrand(i), Trand(i));
end
fprintf('total time: %.0fs\n', toc(tall));

save('figure2_data.mat', 'MS', 'N', 'z', 'b', 'Wex', 'Wrand', 'Wfloor', 'Tex', 'Trand');

%% ------------------------------------------------------------ plotting ----
BLUE = [0 0.447 0.698];  VERM = [0.835 0.369 0];  GREY = [0.54 0.54 0.54];
mm = logspace(log10(MS(1)), log10(MS(end)), 200);

figure(1); clf
set(gcf, 'Color', 'w');
loglog(MS, Wrand, 's-', 'Color', VERM, 'MarkerFaceColor', VERM, ...
       'MarkerSize', 4, 'LineWidth', 1.3); hold on
loglog(MS, Wex, 'o--', 'Color', BLUE, 'MarkerFaceColor', BLUE, ...
       'MarkerSize', 4, 'LineWidth', 1.3);
loglog(mm, 3./mm, '--', 'Color', GREY, 'LineWidth', 2);
yline(Wfloor, ':', 'Color', VERM, 'LineWidth', 2);

xlabel('$m$', 'Interpreter', 'latex')
ylabel('$W_1(\mu_A,\widehat\mu_A)$', 'Interpreter', 'latex')
title('Earth mover''s error of Algorithm 3', 'Interpreter', 'latex')
legend({'random $b$', 'exact $b$ ($\mu_{A,b}=\mu_A$)', '$O(m^{-1})$', ...
        '$W_1(\mu_A,\mu_{A,b})$'}, 'Interpreter', 'latex', 'Location', 'northeast')
xlim([MS(1), MS(end)])

exportgraphics(gcf, 'figure2.png', 'Resolution', 300);
exportgraphics(gcf, 'figure2.png', 'Resolution', 300);