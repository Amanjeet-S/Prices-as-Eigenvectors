function [z, J] = excessDemand(econ, p)
%EXCESSDEMAND Aggregate excess demand and its analytic Jacobian.
%   [z, J] = EXCESSDEMAND(econ, p) returns z(p) = sum_i x_i(p) - supply and
%   J(j,l) = dz_j/dp_l, using
%       dx_ij/dp_l = x_ij*(e_il/income_i - (1-sigma) a_il p_l^(-sigma)/D_i)
%                    - sigma*x_ij/p_j*[j == l].

arguments
    econ (1,1) struct
    p (:,1) double {mustBePositive}
end

[X, income] = demand(econ, p);
z = sum(X, 1)' - econ.supply;

if nargout > 1
    s = econ.sigma;
    D = econ.A * p.^(1 - s);
    dLogX = econ.E ./ income - (1 - s) * econ.A .* (p.^(-s))' ./ D;
    J = X' * dLogX - s * diag(sum(X, 1)' ./ p);
end
end
