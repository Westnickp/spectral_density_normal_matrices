function d = w1_pot(zA, a, zB, b)
%W1_POT  Exact Wasserstein-1 distance between two atomic measures in the
% plane, via the network simplex solver in POT:
%   R. Flamary et al., POT: Python Optimal Transport, J. Mach. Learn. Res.
%   22(78) (2021), 1-8.
%
%   d = w1_pot(zA, a, zB, b)
%
% zA, zB: atom locations (complex vectors); a, b: nonnegative masses
% (normalized here).  Returns min sum_ij F_ij |zA_i - zB_j| over couplings F.

    a = a(:) / sum(a);  b = b(:) / sum(b);
    N = numel(a);  M = numel(b);
    C = abs(zA(:) - zB(:).');                 % N-by-M distance matrix

    aP = to_numpy(a,   N);
    bP = to_numpy(b,   M);
    CP = to_numpy(C.', [N M]);                % C.' so C is flattened by rows

    d = double(py.ot.emd2(aP, bP, CP, pyargs('numItermax', int64(1e8))));
end

function X = to_numpy(x, shp)
%TO_NUMPY  Send x to Python as a flat float64 array, then reshape to shp
% (row-major).  ndmin=1 makes a scalar a length-1 array.
    flat = py.numpy.array(reshape(x, 1, []), ...
                          pyargs('dtype', 'float64', 'ndmin', int64(1)));
    X = py.numpy.reshape(flat, py.tuple(num2cell(int64(shp))));
end