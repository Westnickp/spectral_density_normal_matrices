function rho = jackson_coeffs(m)
%JACKSON_COEFFS  Jackson damping coefficients.
%
%   rho = jackson_coeffs(m)
%
% Returns the (m+1)-by-1 vector rho(j+1) = rho_j, j = 0..m, with
%
%   rho_j = [ (m+2-j) cos(j*pi/(m+2)) + sin(j*pi/(m+2)) cot(pi/(m+2)) ] / (m+2).

    j     = (0:m).';
    alpha = pi / (m+2);
    rho   = ((m+2-j).*cos(j*alpha) + sin(j*alpha)*cot(alpha)) / (m+2);
end
