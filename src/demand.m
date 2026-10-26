function [X, income] = demand(econ, p)
%DEMAND CES demands of every consumer at prices p.
%   [X, income] = DEMAND(econ, p) returns one row of X per consumer, from
%   the Lagrange conditions:
%       x_ij = a_ij p_j^(-sigma) income_i / sum_k a_ik p_k^(1-sigma).

arguments
    econ (1,1) struct
    p (:,1) double {mustBePositive}
end

s = econ.sigma;
income = econ.E * p;
D = econ.A * p.^(1 - s);
X = econ.A .* (p.^(-s))' .* (income ./ D);
end
