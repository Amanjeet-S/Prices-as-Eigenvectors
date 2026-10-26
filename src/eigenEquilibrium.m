function eqm = eigenEquilibrium(econ)
%EIGENEQUILIBRIUM Cobb-Douglas equilibrium as an eigenvector.
%   eqm = EIGENEQUILIBRIUM(econ) solves v = B*v, where v_j = p_j*supply_j
%   and B = A' * E * diag(1./supply) has columns summing to one.
%   Fields: B, lambda (eigenvalues by decreasing modulus), nullity
%   (dim null(B - I)), v (values summing to one) and p (prices, p(1) = 1).
%   v and p are NaN when nullity > 1, as prices are then not unique.

arguments
    econ (1,1) struct
end

if econ.sigma ~= 1
    error('eigenEquilibrium:notCobbDouglas', ...
        'The eigenvector solution needs Cobb-Douglas preferences (sigma = 1).');
end

eqm.B = econ.A' * econ.E ./ econ.supply';
lambda = eig(eqm.B);
[~, order] = sort(abs(lambda), 'descend');
eqm.lambda = lambda(order);

K = null(eqm.B - eye(econ.n));
eqm.nullity = size(K, 2);
if eqm.nullity ~= 1
    eqm.v = NaN(econ.n, 1);
    eqm.p = NaN(econ.n, 1);
    return
end

% Scaling to sum to one fixes the arbitrary sign and scale from null.
eqm.v = K / sum(K);
p = eqm.v ./ econ.supply;
eqm.p = p / p(1);
end
