function [p, history] = newtonEquilibrium(econ, p0, opts)
%NEWTONEQUILIBRIUM Walrasian equilibrium by Newton's method.
%   [p, history] = NEWTONEQUILIBRIUM(econ, p0) solves z(p) = 0 from p0 and
%   returns prices with p(1) = 1. Good 1 is the numeraire, so only
%   z_2..z_n are solved (Walras' law clears market 1). The unknowns are log
%   prices, so dz_i/dtheta_j = J(i,j)*p_j. history records the residual
%   and step norms; iteration stops when the step is below opts.Tol.

arguments
    econ (1,1) struct
    p0 (:,1) double {mustBePositive}
    opts.Tol (1,1) double {mustBePositive} = 1e-13
    opts.MaxIter (1,1) double {mustBeInteger, mustBePositive} = 50
end

theta = log(p0(2:end) / p0(1));
residual = zeros(opts.MaxIter, 1);
stepNorm = zeros(opts.MaxIter, 1);

for k = 1:opts.MaxIter
    p = [1; exp(theta)];
    [z, J] = excessDemand(econ, p);
    f = z(2:end);
    Df = J(2:end, 2:end) .* p(2:end)';
    step = -(Df \ f);
    theta = theta + step;

    residual(k) = norm(f);
    stepNorm(k) = norm(step);
    if stepNorm(k) < opts.Tol
        p = [1; exp(theta)];
        history = table((1:k)', residual(1:k), stepNorm(1:k), ...
            'VariableNames', {'Iteration', 'Residual', 'Step'});
        return
    end
end

error('newtonEquilibrium:noConvergence', ...
    'Newton''s method did not converge in %d iterations.', opts.MaxIter);
end
