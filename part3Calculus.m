%% Part III: Derivatives, implicit functions and stability (course modules 9-13)
% CES equilibrium by Newton's method, Lagrange and implicit-function
% comparative statics, an economy with three equilibria, and Scarf's cycles.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));
outDir = fullfile(root, 'results');
if ~isfolder(outDir)
    mkdir(outDir);
end
fprintf('\n===== Part III: derivatives, implicit functions and stability =====\n');

[shares, endowments, goods] = exampleEconomy();
sigma = 0.5;                                   % goods are complements
econ = makeEconomy(shares, endowments, sigma);

%% 1. The Jacobian and why it is singular (modules 3 and 10)
% Homogeneity gives J*p = 0 and Walras' law gives p'*J = -z', so at an
% equilibrium rank J <= n - 1 and one price is fixed as numeraire.
p0 = ones(econ.n, 1);
[z0, J0] = excessDemand(econ, p0);
Jnum = numericalJacobian(@(p) excessDemand(econ, p), p0);
fprintf('Analytic vs central-difference Jacobian: relative difference %.1e\n', ...
    norm(J0 - Jnum) / norm(J0));
fprintf('||J p|| = %.1e and ||p''J + z''|| = %.1e at p = (1, 1, 1)\n', ...
    norm(J0 * p0), norm(p0' * J0 + z0'));

%% 2. Newton's method in log prices (modules 10 and 13)
% Taylor's theorem predicts step(k+1) ~ C * step(k)^2.
[pStar, history] = newtonEquilibrium(econ, p0);
history.Ratio = [NaN; history.Step(2:end) ./ history.Step(1:end-1).^2];
disp(history);
writetable(history, fullfile(outDir, 'part3_newton.csv'));

[zStar, JStar] = excessDemand(econ, pStar);
[U, Sv, W] = svd(JStar);
fprintf('Equilibrium with sigma = %.1f: %s, max |z| = %.1e\n', ...
    sigma, mat2str(pStar', 6), max(abs(zStar)));
fprintf('Singular values of J(p*): %s\n', mat2str(diag(Sv)', 4));
fprintf('Angle between p* and the right / left null vectors: %.1e / %.1e rad\n', ...
    subspace(W(:, end), pStar), subspace(U(:, end), pStar));

%% 3. Prices move continuously with the elasticity (modules 7 and 8)
% Trace the equilibrium from sigma = 0.25 to 4, starting from the Part I
% eigenvector at sigma = 1 and warm-starting each solve.
sigmas = logspace(log10(0.25), log10(4), 41);
mid = (numel(sigmas) + 1) / 2;                 % sigmas(mid) = 1
cobbDouglas = eigenEquilibrium(makeEconomy(shares, endowments));
pathPrices = zeros(econ.n, numel(sigmas));
pathPrices(:, mid) = cobbDouglas.p;
for k = [mid+1:numel(sigmas), mid-1:-1:1]
    previous = pathPrices(:, k - sign(k - mid));
    pathPrices(:, k) = newtonEquilibrium(makeEconomy(shares, endowments, sigmas(k)), previous);
end
pNearOne = newtonEquilibrium(makeEconomy(shares, endowments, 1 + 1e-6), p0);
fprintf('sigma = 1 + 1e-6 against the eigenvector solution: max difference %.1e\n', ...
    max(abs(pNearOne - cobbDouglas.p)));

%% 4. Lagrange multipliers and the envelope theorem (module 12)
% At demand, grad u = lambda * p with lambda = 1/P, the marginal utility of income.
XStar = demand(econ, pStar);
for i = 1:econ.m
    [~, g] = utility(shares(i, :), XStar(i, :), sigma);
    multipliers = g ./ pStar';
    priceIndex = (shares(i, :) * pStar.^(1 - sigma))^(1 / (1 - sigma));
    fprintf('Household %d: grad u ./ p varies by %.1e; lambda = %.6f and 1/P = %.6f\n', ...
        i, max(multipliers) - min(multipliers), mean(multipliers), 1 / priceIndex);
end

%% 5. Comparative statics with the implicit function theorem (module 11)
% Move tau units of food from household 1 to household 2:
% dp/dtau = -(D_p z)^(-1) * D_tau z, with D_tau z exact since demand is linear in income.
transfer = zeros(size(endowments));
transfer(1, 1) = -1;
transfer(2, 1) = 1;
[X0, income0] = demand(econ, pStar);
dzdtau = (X0 ./ income0)' * (transfer * pStar);
dpIFT = -JStar(2:end, 2:end) \ dzdtau(2:end);
% Household 2 owns no food, so check with a second-order one-sided difference.
dTau = 1e-5;
pOne = newtonEquilibrium(makeEconomy(shares, endowments + dTau * transfer, sigma), pStar);
pTwo = newtonEquilibrium(makeEconomy(shares, endowments + 2 * dTau * transfer, sigma), pStar);
dpResolve = (-3 * pStar(2:end) + 4 * pOne(2:end) - pTwo(2:end)) / (2 * dTau);
fprintf('Food transfer, d(clothing, fuel prices)/dtau: %s (IFT), %s (re-solving)\n', ...
    mat2str(dpIFT', 6), mat2str(dpResolve', 6));
fprintf('Household 2 spends most of its extra income on clothing, so clothing gets dearer.\n');

%% 6. Equilibrium lies on the contract curve (modules 2 and 12)
% Both gradients are parallel to p, so [grad u1; grad u2] is singular:
% the tangency behind the First Welfare Theorem.
edgeShares = [0.6 0.4; 0.3 0.7];
edgeEndowments = [3 1; 1 3];
edge = makeEconomy(edgeShares, edgeEndowments);
edgeEqm = eigenEquilibrium(edge);
XEdge = demand(edge, edgeEqm.p);
[u1, g1] = utility(edgeShares(1, :), XEdge(1, :));
[u2, g2] = utility(edgeShares(2, :), XEdge(2, :));
fprintf('Edgeworth box: |sin(angle between gradients)| = %.1e\n', ...
    abs(det([g1; g2])) / (norm(g1) * norm(g2)));

%% 7. Three equilibria and where the implicit function theorem fails (modules 9 and 11)
% The equilibrium set is the line q = 1 plus the curve
% r(q) = (q^(8/9) - q^(1/9)) / (q - 1). They cross at (7/9, 1), where the
% derivative of z2 vanishes and the set is not a manifold.
r0 = 2^(8/9) - 2^(1/9);
qEq = [0.5; 1; 2];
[zEq, dzEq] = textbookExcessDemand(qEq, r0);
fprintf('Textbook economy, r = %.6f: z2 = %s at q = 1/2, 1, 2\n', r0, mat2str(zEq', 3));
fprintf('dz2/dq = %s, so tatonnement-stable: %s; index sum = %d\n', ...
    mat2str(dzEq', 4), mat2str(dzEq' < 0), sum(sign(-dzEq)));
[~, dzCrossing] = textbookExcessDemand(1, 7/9);
fprintf('At the crossing (7/9, 1): dz2/dq = %.1e and dz2/dr = 1/q - 1 = 0\n', dzCrossing);

qCurve = logspace(log10(0.05), log10(20), 801)';
qCurve(abs(qCurve - 1) < 1e-9) = [];
rCurve = (qCurve.^(8/9) - qCurve.^(1/9)) ./ (qCurve - 1);
[~, dzCurve] = textbookExcessDemand(qCurve, rCurve);

%% 8. Scarf's economy: prices that never settle (modules 2, 4, 11 and 13)
% Tatonnement dp/dt = z(p) conserves ||p||^2 and p1*p2*p3, so paths circle
% the equilibrium (1, 1, 1) and never reach it.
odeOpts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
pStart = [1.2; 1.0; 0.8];
tEnd = 60;
[tScarf, PScarf] = ode45(@(t, p) scarfExcessDemand(p), [0 tEnd], pStart, odeOpts);
fprintf('Scarf: largest change in ||p||^2 is %.1e and in p1*p2*p3 is %.1e over t in [0, %d]\n', ...
    max(abs(sum(PScarf.^2, 2) - sum(pStart.^2))), max(abs(prod(PScarf, 2) - prod(pStart))), tEnd);

% The Jacobian at (1, 1, 1) is skew-symmetric: eigenvalues +-i*sqrt(3)/4, period 8*pi/sqrt(3).
JScarf = numericalJacobian(@scarfExcessDemand, ones(3, 1));
muScarf = eig(JScarf);
fprintf('||J + J''|| = %.1e; eigenvalues have real parts below %.1e and imaginary parts +-%.4f (theory %.4f)\n', ...
    norm(JScarf + JScarf'), max(abs(real(muScarf))), max(imag(muScarf)), sqrt(3)/4);
fprintf('Period of a small cycle: %.4f (theory 8*pi/sqrt(3) = %.4f)\n', ...
    scarfCyclePeriod([1.001; 1; 0.999]), 8 * pi / sqrt(3));

% The Part I economy has gross substitutes, so tatonnement converges
% (Arrow, Block and Hurwicz, 1959). Distances use the equilibrium ray at the same norm.
cdEcon = makeEconomy(shares, endowments);
[tCD, PCD] = ode45(@(t, p) excessDemand(cdEcon, p), [0 tEnd], pStart, odeOpts);
cdTarget = cobbDouglas.p / norm(cobbDouglas.p) * norm(pStart);
scarfTarget = ones(3, 1) / sqrt(3) * norm(pStart);
distanceCD = vecnorm(PCD' - cdTarget);
distanceScarf = vecnorm(PScarf' - scarfTarget);
fprintf('Distance to equilibrium at t = %d: %.1e (Part I economy), %.3f (Scarf)\n', ...
    tEnd, distanceCD(end), distanceScarf(end));

%% Figure
fig = figure('Name', 'Part III: calculus', 'Position', [50 50 1400 820]);
tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
blue = [0 0.447 0.741];
red = [0.850 0.325 0.098];

% (a) Newton's method
nexttile;
semilogy(history.Iteration, history.Residual, 'o-', history.Iteration, history.Step, 's-', ...
    'LineWidth', 1.2);
xlabel('Iteration');
ylabel('Norm');
legend({'||z(p_k)||', '||step||'}, 'Location', 'southwest');
title('(a) Newton converges quadratically');

% (b) Equilibrium prices against the elasticity of substitution
nexttile;
loglog(sigmas, pathPrices(2:3, :)', 'LineWidth', 1.5);
hold on
loglog(1, cobbDouglas.p(2), 'ko', 1, cobbDouglas.p(3), 'ko', 'MarkerFaceColor', 'k');
hold off
xlabel('Elasticity of substitution \sigma');
ylabel('Price (food = 1)');
legend([goods(2:3); "Eigenvector solution (\sigma = 1)"], 'Location', 'northeast');
title('(b) Prices vary continuously with \sigma');

% (c) Edgeworth box
nexttile;
supply = edge.supply;
[G1, G2] = meshgrid(linspace(0.02, supply(1) - 0.02, 200), linspace(0.02, supply(2) - 0.02, 200));
U1 = reshape(utility(edgeShares(1, :), [G1(:), G2(:)]), size(G1));
U2 = reshape(utility(edgeShares(2, :), [supply(1) - G1(:), supply(2) - G2(:)]), size(G1));
contour(G1, G2, U1, [u1 u1], 'LineColor', blue, 'LineWidth', 1.2);
hold on
contour(G1, G2, U2, [u2 u2], 'LineColor', red, 'LineWidth', 1.2);
ratio1 = edgeShares(1, 1) / edgeShares(1, 2);
ratio2 = edgeShares(2, 1) / edgeShares(2, 2);
x11 = linspace(0, supply(1), 200);
contractCurve = ratio2 * supply(2) * x11 ./ (ratio1 * (supply(1) - x11) + ratio2 * x11);
plot(x11, contractCurve, 'k-', 'LineWidth', 1);
budgetX = linspace(0, supply(1), 2);
budgetY = edgeEndowments(1, 2) - edgeEqm.p(1) / edgeEqm.p(2) * (budgetX - edgeEndowments(1, 1));
plot(budgetX, budgetY, '--', 'Color', [0.4 0.4 0.4]);
plot(edgeEndowments(1, 1), edgeEndowments(1, 2), 'ks', 'MarkerFaceColor', 'w');
plot(XEdge(1, 1), XEdge(1, 2), 'ko', 'MarkerFaceColor', 'k');
hold off
axis([0 supply(1) 0 supply(2)]);
axis square
xlabel('Good 1 held by household 1');
ylabel('Good 2 held by household 1');
legend({'Household 1 indifference', 'Household 2 indifference', 'Contract curve', ...
    'Budget line', 'Endowment', 'Equilibrium'}, 'Location', 'southoutside', ...
    'NumColumns', 2, 'FontSize', 7);
title('(c) Equilibrium lies on the contract curve');

% (d) The equilibrium set of the textbook economy
nexttile;
stableCurve = rCurve;
stableCurve(dzCurve >= 0) = NaN;
unstableCurve = rCurve;
unstableCurve(dzCurve < 0) = NaN;
semilogy(stableCurve, qCurve, '-', 'Color', blue, 'LineWidth', 1.8);
hold on
semilogy(unstableCurve, qCurve, '-', 'Color', red, 'LineWidth', 1.8);
semilogy([7/9 0.85], [1 1], '-', 'Color', blue, 'LineWidth', 1.8);
semilogy([0.70 7/9], [1 1], '-', 'Color', red, 'LineWidth', 1.8);
xline(r0, ':', 'r = 2^{8/9} - 2^{1/9}', 'LabelVerticalAlignment', 'top');
semilogy(r0 * [1 1 1], qEq, 'ko', 'MarkerFaceColor', 'k');
semilogy(7/9, 1, 'kp', 'MarkerSize', 9, 'MarkerFaceColor', 'y');
text(7/9, 1.3, ' r = 7/9', 'FontSize', 8);
hold off
xlim([0.70 0.85]);
ylim([0.05 20]);
xlabel('Endowment parameter r');
ylabel('Equilibrium price ratio q = p_2/p_1');
legend({'Stable equilibria', 'Unstable equilibria'}, 'Location', 'northeast');
title('(d) Three equilibria; the IFT fails at the star');

% (e) A Scarf cycle on its sphere
nexttile;
radius = norm(pStart);
[Xs, Ys, Zs] = sphere(60);
surf(radius * Xs, radius * Ys, radius * Zs, 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
    'FaceColor', [0.6 0.6 0.6]);
hold on
plot3(PScarf(:, 1), PScarf(:, 2), PScarf(:, 3), 'Color', red, 'LineWidth', 1.5);
plot3(scarfTarget(1), scarfTarget(2), scarfTarget(3), 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
hold off
axis equal
axis([0.6 1.4 0.6 1.4 0.6 1.4]);
view(135, 25);
xlabel('p_1');
ylabel('p_2');
zlabel('p_3');
title('(e) Scarf: prices circle the equilibrium');

% (f) Distance to equilibrium under tatonnement
nexttile;
semilogy(tCD, distanceCD, 'LineWidth', 1.5, 'Color', blue);
hold on
semilogy(tScarf, distanceScarf, 'LineWidth', 1.5, 'Color', red);
hold off
xlim([0 20]);
ylim([1e-14 10]);
xlabel('Time t');
ylabel('Distance to equilibrium');
legend({'Part I economy (gross substitutes)', 'Scarf economy (complements)'}, 'Location', 'east');
title('(f) Tatonnement need not converge');

exportgraphics(fig, fullfile(outDir, 'part3_calculus.png'), 'Resolution', 200);
