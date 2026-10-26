function econ = makeEconomy(shares, endowments, sigma)
%MAKEECONOMY Pure exchange economy with CES consumers.
%   econ = MAKEECONOMY(shares, endowments, sigma) builds an economy with
%   one row per consumer: preference weights (rows sum to one) and
%   endowments. sigma is the elasticity of substitution; sigma = 1 (the
%   default) is Cobb-Douglas.

arguments
    shares (:,:) double {mustBeNonnegative}
    endowments (:,:) double {mustBeNonnegative}
    sigma (1,1) double {mustBePositive} = 1
end

if ~isequal(size(shares), size(endowments))
    error('makeEconomy:size', 'Shares and endowments must have the same size.');
end
if any(abs(sum(shares, 2) - 1) > 1e-12)
    error('makeEconomy:shares', 'Each row of shares must sum to one.');
end
if any(sum(endowments, 1) <= 0)
    error('makeEconomy:supply', 'Every good must be owned by someone.');
end
if any(sum(endowments, 2) <= 0)
    error('makeEconomy:income', 'Every consumer must own something.');
end

econ.A = shares;
econ.E = endowments;
econ.sigma = sigma;
econ.supply = sum(endowments, 1)';   % aggregate endowment, n-by-1
[econ.m, econ.n] = size(shares);
end
