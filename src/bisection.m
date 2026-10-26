function [root, history] = bisection(f, lo, hi, tol)
%BISECTION Zero of a continuous function by repeated halving.
%   [root, history] = BISECTION(f, lo, hi) needs f(lo) and f(hi) of
%   opposite sign. The nested intervals shrink to a zero, which is the
%   intermediate value theorem made constructive. Each row of history is
%   one interval [lo, hi].

arguments
    f function_handle
    lo (1,1) double
    hi (1,1) double
    tol (1,1) double {mustBePositive} = 4 * eps
end

flo = f(lo);
if sign(flo) * sign(f(hi)) >= 0
    error('bisection:noSignChange', 'f must change sign on [lo, hi].');
end

history = [lo, hi];
while hi - lo > tol * max(1, abs(lo))
    mid = lo + (hi - lo) / 2;
    fmid = f(mid);
    if fmid == 0
        lo = mid;
        hi = mid;
    elseif sign(fmid) == sign(flo)
        lo = mid;
        flo = fmid;
    else
        hi = mid;
    end
    history(end + 1, :) = [lo, hi]; %#ok<AGROW>
end
root = lo + (hi - lo) / 2;
end
