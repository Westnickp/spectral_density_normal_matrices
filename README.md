# Supplementary code for "Spectral Density Estimation for Normal Matrices"

This folder contains MATLAB implementations of Algorithms 1–3 of the paper "Spectral density estimation for normal matrices" (https://arxiv.org/abs/2605.31430) and two scripts that reproduce Figures 1 and 2.

## Requirements

- MATLAB R2025b.
- `figure1.m` and the algorithm files run entirely in MATLAB.
- `figure2.m` requires Python with the package POT (Python Optimal Transport), which computes the Wasserstein distances in Figure 2:

  > R. Flamary et al., POT: Python Optimal Transport, *J. Mach. Learn. Res.* 22(78) (2021), 1–8.

  One way to ensure this is installed:

  1. In MATLAB, run `pe = pyenv` and note `pe.Executable`, the Python interpreter MATLAB uses.
  2. In a terminal, install POT into that interpreter:
     ```bash
     /path/to/that/python3 -m pip install pot
     ```
  3. Restart MATLAB and check that
     ```matlab
     py.importlib.import_module('ot')
     ```
     runs without error.

## Quick start

Place all files in one folder, make it the current folder in MATLAB, and run

```matlab
figure1      % writes figure1.png
figure2      % writes figure2.png and figure2_data.mat
```

Each script prints its timings to the command window. Both scripts fix the random seed (`rng`). `figure2` may require installation of Python Optimal Transport (POT); see [Requirements](#requirements).

## Files

### Scripts

| File | Description |
|---|---|
| `figure1.m` | Reproduces Figure 1 (density estimates for the disk/quadrant matrix and the "Anscombosaurus" matrix, m = 16, 32, 64, 128). |
| `figure2.m` | Reproduces Figure 2 (earth mover's error of Algorithm 3 versus m = 2, 4, ..., 128, with exact and random starting vectors). |

### Algorithms

| File | Description |
|---|---|
| `cheb_moments.m` | Algorithm 1 (ChebMoments): mixed normalized Chebyshev moments `Gamma` of `mu_{A,b}`, using 4m products with `A`, `A'`. |
| `jackson_coeffs.m` | Jackson damping coefficients `rho_j`, eq. (jacksoncoeffs). |
| `sde.m` | Algorithms 2 and 3 (WSpecDens / SpecDens): returns the damped moments `M(j+1,k+1) = rho_j rho_k Gamma_jk` that define the estimate. With no `b` supplied, `b` is drawn uniformly from the unit sphere (Algorithm 3). |
| `kpm_density.m` | Evaluates the estimated density `p_KPM(x,y) / (sqrt(1-x^2) sqrt(1-y^2))` on a tensor grid. |
| `sde_atoms.m` | The (m+1)^2-atom Gauss–Chebyshev discretization of the estimate used in the proof of Theorem 1.1. It has the same Chebyshev moments through degree m as the continuous estimate. |
| `cheb_evaluation_1d.m` | Normalized Chebyshev polynomials `Ttilde_0, ..., Ttilde_m`. |

### Utilities

| File | Description |
|---|---|
| `w1_pot.m` | Exact Wasserstein-1 (earth mover's) distance between two atomic measures in the plane, computed by the network simplex solver in POT (`ot.emd2`). Used by `figure2.m`. |
| `preprocess.m` | Estimates a shift `c2` and scale `c1` so that `c1*(A - c2*I)` has norm < 1, using a few products with `A` and `A'` (cf. the norm assumption `norm(A) <= 1` in the paper). Not needed for the figures, whose spectra are scaled exactly. |

### Data

| File | Description |
|---|---|
| `pixels.csv` | 27,077 points (x,y) in [0,1]^2 tracing the "Anscombosaurus"; the eigenvalues of the second matrix in Figure 1. |

## Conventions

- **m is the maximum polynomial degree**, as in the paper. Moment and coefficient matrices are therefore (m+1)-by-(m+1), and MATLAB index `(j+1,k+1)` corresponds to the paper's index `(j,k)`.
- A complex number x + iy is identified with the point (x,y). The estimate lives on the square [-1,1]^2, which requires `norm(A) <= 1`.
- Grid outputs of `kpm_density.m` follow `meshgrid`: rows correspond to y and columns to x, so they can be passed directly to `imagesc` or `contourf`.
- All functions access `A` only through products `A*v` and `A'*v`.

## Example

```matlab
n = 2000;
z = randn(n,1) + 1i*randn(n,1);                    % eigenvalues

c2 = mean(z);  c1 = 0.95/max(abs(z - c2));         % shift and scale so
A  = spdiags(c1*(z - c2), 0, n, n);                % that norm(A) < 1

M = sde(A, 50);                                    % Algorithm 3, m = 50

x = linspace(-0.99, 0.99, 400);
F = kpm_density(M, x, x);                          % estimated density
imagesc(x, x, F), axis xy equal tight, colormap(flipud(gray))
```

The point (x,y) of the plot corresponds to the eigenvalue location `c2 + (x + 1i*y)/c1` of the original matrix.

## Reproducibility notes

- `figure1.m`: `rng(1)` fixes both test spectra and the random starting vectors.
- `figure2.m`: `rng(1)` fixes the test spectrum and `rng(2)` fixes the random `b`, so that `b` is independent of the test matrix.
