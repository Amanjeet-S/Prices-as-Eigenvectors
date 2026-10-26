function z = scarfExcessDemand(p)
%SCARFEXCESSDEMAND Excess demand in Scarf's (1960) three-good economy.
%   Consumer i owns one unit of good i and has Leontief utility
%   min(x_i, x_{i+1}), indices taken mod 3. They spend their income p_i on
%   equal amounts of goods i and i+1, so excess demand for good i is
%       z_i = p_i/(p_i + p_{i+1}) + p_{i-1}/(p_{i-1} + p_i) - 1.
%   p may be a 3-by-1 price vector or a 3-by-k matrix of them.

arguments
    p (3,:) double {mustBePositive}
end

next = circshift(p, -1, 1);   % row i holds p_{i+1}
prev = circshift(p, 1, 1);    % row i holds p_{i-1}
z = p ./ (p + next) + prev ./ (prev + p) - 1;
end
