function T = scarfCyclePeriod(p0)
%SCARFCYCLEPERIOD Length of one tatonnement cycle in Scarf's economy.
%   T = SCARFCYCLEPERIOD(p0) integrates dp/dt = z(p) from p0 with ode45 and
%   returns the first time after t = 1 at which p_1 comes back to p0(1)
%   moving in the same direction as at the start: one full cycle.

arguments
    p0 (3,1) double {mustBePositive}
end

z0 = scarfExcessDemand(p0);
opts = odeset('RelTol', 1e-11, 'AbsTol', 1e-13, ...
    'Events', @(t, p) returnToStart(p, p0(1), sign(z0(1))));
sol = ode45(@(t, p) scarfExcessDemand(p), [0 200], p0, opts);
T = sol.xe(find(sol.xe > 1, 1));
if isempty(T)
    error('scarfCyclePeriod:noReturn', 'The path did not return before t = 200.');
end
end

function [value, isTerminal, direction] = returnToStart(p, p1, dir)
value = p(1) - p1;
isTerminal = 0;
direction = dir;
end
