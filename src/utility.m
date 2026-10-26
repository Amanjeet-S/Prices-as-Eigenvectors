function [u, grad] = utility(a, x, sigma)
%UTILITY CES utility of each row of x, and its gradient.
%   u = ( sum_j a_j^(1/sigma) x_j^((sigma-1)/sigma) )^(sigma/(sigma-1)),
%   with gradient u^(1/sigma) a_j^(1/sigma) x_j^(-1/sigma). sigma = 1 uses
%   the limit prod_j (x_j/a_j)^a_j, so u is continuous in sigma.

arguments
    a (1,:) double {mustBeNonnegative}
    x (:,:) double {mustBePositive}
    sigma (1,1) double {mustBePositive} = 1
end

if sigma == 1
    u = prod((x ./ a).^a, 2);      % a zero weight gives Inf^0 = 1
else
    r = (sigma - 1) / sigma;
    u = sum(a.^(1/sigma) .* x.^r, 2).^(1/r);
end

if nargout > 1
    grad = u.^(1/sigma) .* a.^(1/sigma) .* x.^(-1/sigma);
end
end
