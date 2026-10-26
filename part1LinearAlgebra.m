%% Part I: Prices as eigenvectors (course modules 0-4)
% Cobb-Douglas equilibrium as an eigenvector, a counterexample to
% uniqueness, and the same prices as the limit of a sequence.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));
outDir = fullfile(root, 'results');
if ~isfolder(outDir)
    mkdir(outDir);
end

[shares, endowments, goods] = exampleEconomy();
econ = makeEconomy(shares, endowments);
fprintf('\n===== Part I: prices as eigenvectors =====\n');

%% 1. The expenditure-flow matrix B (modules 1, 3 and 4)
% Market clearing is v = B*v with v_j = p_j * supply_j. The columns of B
% sum to one, so the equilibrium values are an eigenvector for eigenvalue 1.
eqm = eigenEquilibrium(econ);
fprintf('Column sums of B: %s\n', mat2str(sum(eqm.B, 1), 12));
fprintf('Eigenvalues of B: %s\n', mat2str(eqm.lambda.', 6));
disp('Row reduction of B - I leaves a zero row, so rank(B - I) = 2:');
disp(rref(eqm.B - eye(econ.n)));
fprintf('null(B - I) has dimension %d, so relative prices are unique\n', eqm.nullity);

%% 2. Checking the eigenvector against the definition of equilibrium
z = excessDemand(econ, eqm.p);
prices = table(goods, eqm.p, z, 'VariableNames', {'Good', 'Price', 'ExcessDemand'});
disp(prices);
writetable(prices, fullfile(outDir, 'part1_prices.csv'));

% Row reduction never leaves the field Q, so the prices are rational (module 0).
fprintf('As fractions: %s\n', strtrim(rats(eqm.p')));

% Walras' law: z(p) is orthogonal to p at every price vector (module 2).
rng(2026);
pRandom = 0.5 + rand(econ.n, 1);
fprintf('Walras'' law at a random price vector: p''z(p) = %.1e\n', ...
    pRandom' * excessDemand(econ, pRandom));

%% 3. When the hypothesis fails: an economy of two islands
% Two groups that never trade make B block diagonal, so null(B - I) is
% two-dimensional and relative prices across the islands are undetermined.
islandShares = [0.7 0.3 0 0; 0.4 0.6 0 0; 0 0 0.5 0.5; 0 0 0.2 0.8];
islandEndowments = [2 1 0 0; 1 3 0 0; 0 0 1 2; 0 0 3 1];
islands = eigenEquilibrium(makeEconomy(islandShares, islandEndowments));
fprintf('Two islands: eigenvalue 1 appears %d times and null(B - I) has dimension %d\n', ...
    sum(abs(islands.lambda - 1) < 1e-10), islands.nullity);

%% 4. The same prices as the limit of a sequence (modules 4-6)
% Each round of spending maps v to B*v, and
% ||v_k - v*|| <= cond(V) * |lambda_2|^k * ||v_0 - v*||. Stop at 18 rounds,
% before the error reaches rounding level.
[V, ~] = eig(eqm.B);
lambda2 = abs(eqm.lambda(2));
nSteps = 18;
v = ones(econ.n, 1) / econ.n;
gap = zeros(nSteps + 1, 1);
for k = 0:nSteps
    gap(k + 1) = norm(v - eqm.v);
    v = eqm.B * v;
end
bound = cond(V) * lambda2.^(0:nSteps)' * gap(1);
fprintf('|lambda_2| = %.4f and cond(V) = %.3f; the bound holds at every step: %s\n', ...
    lambda2, cond(V), string(all(gap <= bound * (1 + 1e-9) + 1e-15)));
fprintf('After %d rounds of spending the values are within %.1e of equilibrium\n', ...
    nSteps, gap(end));

%% Figure
fig = figure('Name', 'Part I: prices as eigenvectors', 'Position', [100 100 1000 420]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
theta = linspace(0, 2*pi, 400);
plot(cos(theta), sin(theta), 'Color', [0.6 0.6 0.6]);
hold on
plot(real(eqm.lambda), imag(eqm.lambda), 'o', 'MarkerSize', 6, ...
    'MarkerFaceColor', [0 0.447 0.741], 'MarkerEdgeColor', 'k');
text(1.06, 0.12, '\lambda_1 = 1');
text(0.25, -0.18, sprintf('\\lambda_2 = %.3f', eqm.lambda(2)));
text(0.25, 0.18, sprintf('\\lambda_3 = %.3f', eqm.lambda(3)));
hold off
axis equal
xlim([-1.2 1.5]);
ylim([-1.2 1.2]);
xlabel('Re \lambda');
ylabel('Im \lambda');
title('(a) Eigenvalues of B and the unit circle');

nexttile;
semilogy(0:nSteps, gap, 'o-', 'LineWidth', 1.2);
hold on
semilogy(0:nSteps, bound, '--', 'LineWidth', 1.2);
hold off
xlabel('Round of spending k');
ylabel('||v_k - v^*||');
legend({'Distance to equilibrium', 'cond(V) |\lambda_2|^k ||v_0 - v^*||'}, 'Location', 'northeast');
title('(b) The sequence v_k converges at rate |\lambda_2|');

exportgraphics(fig, fullfile(outDir, 'part1_eigenvectors.png'), 'Resolution', 200);
