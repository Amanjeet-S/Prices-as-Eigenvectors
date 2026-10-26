%% Part II: Existence, sequences and floating point (course modules 5-9)
% A two-good economy with prices (t, 1-t): existence of a zero of z_1,
% bisection as nested intervals, and the accuracy of finite differences.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));
outDir = fullfile(root, 'results');
if ~isfolder(outDir)
    mkdir(outDir);
end

edgeShares = [0.6 0.4; 0.3 0.7];
edgeEndowments = [3 1; 1 3];
edge = makeEconomy(edgeShares, edgeEndowments);
z1 = @(t) firstExcessDemand(edge, t);
fprintf('\n===== Part II: existence, sequences and floating point =====\n');

%% 1. Existence from the intermediate value theorem (modules 7 and 8)
% z_1(t) = S/t - (S + R) tends to +inf as t -> 0 and to -R as t -> 1,
% so it has a zero, which is unique because z_1 is decreasing: t* = S/(S + R).
S = edgeShares(:, 1)' * edgeEndowments(:, 2);
R = (1 - edgeShares(:, 1))' * edgeEndowments(:, 1);
tExact = S / (S + R);
fprintf('z_1 near the edges of the simplex: z_1(1e-8) = %.3g and z_1(1 - 1e-12) = %.4f\n', ...
    z1(1e-8), z1(1 - 1e-12));

%% 2. The equilibrium as the limit of nested intervals (modules 5 and 6)
% Interval widths halve each step, so the midpoints form a Cauchy sequence.
[tBisect, intervals] = bisection(z1, 1e-6, 1 - 1e-6);
tFzero = fzero(z1, [1e-6, 1 - 1e-6]);
fprintf('Bisection: t* = %.15f after %d halvings\n', tBisect, size(intervals, 1) - 1);
fprintf('Closed form: t* = %.15f; fzero: %.15f\n', tExact, tFzero);
fprintf('Equilibrium price ratio p1/p2 = %.6f (closed form S/R = %.6f)\n', ...
    tBisect / (1 - tBisect), S / R);

%% 3. A derivative in floating point is not a limit in R (modules 5 and 9)
% Truncation error is O(h) forward and O(h^2) central; rounding error is
% O(eps/h). The best steps are near sqrt(eps) and eps^(1/3).
slopeExact = -S / tExact^2;
h = 10.^-(1:15)';
forwardErr = abs(arrayfun(@(s) (z1(tExact + s) - z1(tExact)) / s, h) - slopeExact);
centralErr = abs(arrayfun(@(s) (z1(tExact + s) - z1(tExact - s)) / (2 * s), h) - slopeExact);
[bestForward, iF] = min(forwardErr);
[bestCentral, iC] = min(centralErr);
fprintf('Best forward step h = %.0e (error %.1e); best central step h = %.0e (error %.1e)\n', ...
    h(iF), bestForward, h(iC), bestCentral);
fprintf('Theory: sqrt(eps) = %.1e and eps^(1/3) = %.1e\n', sqrt(eps), eps^(1/3));

%% Figure
fig = figure('Name', 'Part II: existence and sequences', 'Position', [100 100 1300 400]);
tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
tGrid = linspace(0.1, 0.95, 300);
plot(tGrid, arrayfun(z1, tGrid), 'LineWidth', 1.5);
hold on
yline(0, '-', 'Color', [0.6 0.6 0.6]);
plot(tExact, 0, 'ko', 'MarkerFaceColor', 'k');
hold off
xlabel('Price of good 1 on the simplex, t');
ylabel('Excess demand z_1(t)');
title('(a) A sign change forces a zero');

nexttile;
hold on
for k = 0:12
    barStart = intervals(k + 1, 1);
    barWidth = intervals(k + 1, 2) - barStart;
    rectangle('Position', [barStart, k - 0.3, barWidth, 0.6], ...
        'FaceColor', [0 0.447 0.741], 'EdgeColor', 'none');
end
xline(tExact, ':', 'Color', [0.3 0.3 0.3]);
hold off
set(gca, 'YDir', 'reverse');
xlim([0 1]);
ylim([-0.5 12.5]);
xlabel('t');
ylabel('Bisection step k');
title('(b) Nested intervals shrink to t^*');
box on

nexttile;
loglog(h, forwardErr, 'o-', h, centralErr, 's-', 'LineWidth', 1.2);
hold on
loglog(h, forwardErr(1) * h / h(1), ':', h, centralErr(1) * (h / h(1)).^2, ':', ...
    h, 4 * eps ./ h, ':', 'Color', [0.5 0.5 0.5]);
hold off
ylim([1e-12 1e2]);
xlabel('Step h');
ylabel('|error in dz_1/dt|');
legend({'Forward difference', 'Central difference', 'O(h), O(h^2), O(eps/h)'}, ...
    'Location', 'north');
title('(c) Truncation against rounding');

exportgraphics(fig, fullfile(outDir, 'part2_existence.png'), 'Resolution', 200);

%% Local functions
function out = firstExcessDemand(econ, t)
% Excess demand for good 1 at the prices (t, 1 - t).
z = excessDemand(econ, [t; 1 - t]);
out = z(1);
end
