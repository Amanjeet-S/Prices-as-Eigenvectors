function J = numericalJacobian(f, p)
%NUMERICALJACOBIAN Central-difference Jacobian of f at p.
%   The step h = eps^(1/3)*max(1, |p_j|) balances O(h^2) truncation error
%   against O(eps/h) rounding error.

arguments
    f function_handle
    p (:,1) double
end

n = numel(p);
J = zeros(numel(f(p)), n);
for j = 1:n
    h = eps^(1/3) * max(1, abs(p(j)));
    step = zeros(n, 1);
    step(j) = h;
    J(:, j) = (f(p + step) - f(p - step)) / (2 * h);
end
end
