%FIGURE1  Reproduces Figure 1: spectral density estimation for two matrices.
%
% Row 1: A_1 = diag(z), n = 150,000.  100,000 eigenvalues uniform in the disk
%        of radius R = 0.9, plus 50,000 uniform in the part of that disk lying
%        in the second and fourth quadrants.
% Row 2: A_2 = diag(z), n = 27,077, eigenvalues tracing the "Anscombosaurus"
%        (pixels.csv), centered at their mean and scaled to radius 0.95.
%
% Column 1 shows the eigenvalues and columns 2-5 show the density output by
% Algorithm 3 (sde.m) for m = 16, 32, 64, 128.  For each matrix a single
% random starting vector b is drawn and used for all four values of m.
% All panels show the square [-1,1] + i[-1,1].
%
% Requires: sde.m, cheb_moments.m, jackson_coeffs.m, cheb_evaluation_1d.m,
%           kpm_density.m, pixels.csv.
% Output:   figure1.png

clear; close all
rng(1);                                  % fixes the eigenvalues and b

MS    = [16 32 64 128];                  % maximum polynomial degrees
NGRID = 500;                             % evaluation grid per axis
xg    = linspace(-1, 1, NGRID+1);
xg    = (xg(1:end-1) + xg(2:end)) / 2;   % cell centres; avoids x = +-1

%% ------------------------------------------------------- test matrices ----
% A_1: disk + second/fourth quadrants
R  = 0.9;
n0 = 100000;  nq = 50000;
r  = R*sqrt(rand(n0,1));  th = 2*pi*rand(n0,1);
z1 = r.*exp(1i*th);
r  = R*sqrt(rand(nq,1));
th = pi/2 + pi*floor(2*rand(nq,1)) + (pi/2)*rand(nq,1);   % in QII or QIV
z1 = [z1; r.*exp(1i*th)];

% A_2: Anscombosaurus
D  = readmatrix('pixels.csv');           % columns x, y in [0,1]
z2 = conj(D(:,1) + 1i*D(:,2));
z2 = z2 - mean(z2);
z2 = 0.95 * z2 / max(abs(z2));

zs = {z1, z2};

%% -------------------------------------------------------------- figure ----
figure(1); clf
set(gcf, 'Units', 'inches', 'Position', [1 1 11 4.4], 'Color', 'w');
tl = tiledlayout(2, 1+numel(MS), 'TileSpacing', 'compact', 'Padding', 'compact');
cmap = flipud(gray(256));
ttotal = tic;

for row = 1:2
    z = zs{row};
    n = numel(z);
    A = spdiags(z, 0, n, n);
    b = randn(n,1) + 1i*randn(n,1);      % one Haar-random b per matrix

    % eigenvalues, drawn as a scatter plot with transparency
    ax = nexttile;
    scatter(ax, real(z), imag(z), 0.5, 'k', 'filled', ...
            'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', 'none');
    format_axes(ax)

    % Algorithm 3 for each m
    for m = MS
        t = tic;
        M = sde(A, m, b);
        F = kpm_density(M, xg, xg);
        fprintf('row %d  n = %6d  m = %3d  time = %.2fs\n', row, n, m, toc(t));

        ax = nexttile;
        imagesc(xg, xg, max(F, 0));
        colormap(ax, cmap);
        format_axes(ax)
    end
end
fprintf('total time (including plotting): %.2fs\n', toc(ttotal));

exportgraphics(gcf, 'figure1.png', 'Resolution', 300);
exportgraphics(gcf, 'figure1.png', 'Resolution', 300);

%% ------------------------------------------------------------- helpers ----
function format_axes(ax)
    axis(ax, 'xy');  axis(ax, 'equal');  axis(ax, [-1 1 -1 1]);
    set(ax, 'XTick', -1:0.5:1, 'YTick', -1:0.5:1, ...
            'XTickLabel', [], 'YTickLabel', [], 'Box', 'on', 'Layer', 'top');
end
