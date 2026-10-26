function [z2, dz2dq] = textbookExcessDemand(q, r)
%TEXTBOOKEXCESSDEMAND A two-good economy with three equilibria.
%   Mas-Colell, Whinston and Green (1995), Section 15.B. With q = p_2/p_1,
%       z2(q, r) = q^(-1/9) - q^(-8/9) + r*(1/q - 1),
%   and dz2dq is its derivative in q. For r = 2^(8/9) - 2^(1/9) the
%   equilibria are q = 1/2, 1 and 2.

arguments
    q double {mustBePositive}
    r double
end

z2 = q.^(-1/9) - q.^(-8/9) + r .* (1 ./ q - 1);
dz2dq = -q.^(-10/9) / 9 + 8 * q.^(-17/9) / 9 - r ./ q.^2;
end
