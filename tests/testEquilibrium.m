function tests = testEquilibrium
%TESTEQUILIBRIUM Every result the project relies on, checked numerically.
%   Run from the project folder with:  results = runtests('tests');

tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'src')));
[shares, endowments] = exampleEconomy();
testCase.TestData.shares = shares;
testCase.TestData.endowments = endowments;
end

%% Linear algebra

function testWalrasLawAndHomogeneity(testCase)
p = [0.7; 1.9; 1.2];
for sigma = [0.5, 1, 2]
    econ = randomEconomy(1, 4, 3, sigma);
    z = excessDemand(econ, p);
    verifyEqual(testCase, p' * z, 0, 'AbsTol', 1e-12);
    verifyEqual(testCase, excessDemand(econ, 3 * p), z, 'AbsTol', 1e-12);
end
end

function testEigenvectorClearsMarkets(testCase)
econ = makeEconomy(testCase.TestData.shares, testCase.TestData.endowments);
eqm = eigenEquilibrium(econ);
verifyEqual(testCase, sum(eqm.B, 1), ones(1, 3), 'AbsTol', 1e-14);
verifyEqual(testCase, eqm.lambda(1), 1, 'AbsTol', 1e-12);
verifyEqual(testCase, eqm.nullity, 1);
verifyEqual(testCase, excessDemand(econ, eqm.p), zeros(3, 1), 'AbsTol', 1e-12);
% Row reduction over the rationals gives exact prices (1, 507/329, 174/329).
verifyEqual(testCase, eqm.p, [1; 507/329; 174/329], 'RelTol', 1e-12);
end

function testReducibleEconomyHasNoUniquePrices(testCase)
shares = [0.7 0.3 0 0; 0.4 0.6 0 0; 0 0 0.5 0.5; 0 0 0.2 0.8];
endowments = [2 1 0 0; 1 3 0 0; 0 0 1 2; 0 0 3 1];
eqm = eigenEquilibrium(makeEconomy(shares, endowments));
verifyEqual(testCase, eqm.nullity, 2);
verifyTrue(testCase, all(isnan(eqm.p)));
end

function testPowerIterationBound(testCase)
eqm = eigenEquilibrium(makeEconomy(testCase.TestData.shares, testCase.TestData.endowments));
[V, ~] = eig(eqm.B);
v = ones(3, 1) / 3;
gap0 = norm(v - eqm.v);
for k = 0:18                     % beyond this the error sits at rounding level

    verifyLessThanOrEqual(testCase, norm(v - eqm.v), ...
        cond(V) * abs(eqm.lambda(2))^k * gap0 + 1e-15);
    v = eqm.B * v;
end
end

%% Real analysis

function testBisectionMatchesClosedForm(testCase)
edge = makeEconomy([0.6 0.4; 0.3 0.7], [3 1; 1 3]);
z1 = @(t) [1 0] * excessDemand(edge, [t; 1 - t]);
[t, intervals] = bisection(z1, 1e-6, 1 - 1e-6);
verifyEqual(testCase, t, 1.5 / 3.4, 'AbsTol', 1e-14);
% The intervals are nested and each is half as wide as the one before.
verifyTrue(testCase, all(intervals(2:end, 1) >= intervals(1:end-1, 1)));
verifyTrue(testCase, all(intervals(2:end, 2) <= intervals(1:end-1, 2)));
widths = diff(intervals, 1, 2);
verifyLessThanOrEqual(testCase, widths(2:end), widths(1:end-1) / 2 + eps);
end

function testCesDemandTendsToCobbDouglas(testCase)
shares = testCase.TestData.shares;
endowments = testCase.TestData.endowments;
p = [1; 2; 0.5];
XCobbDouglas = demand(makeEconomy(shares, endowments), p);
for delta = [1e-3, 1e-5, 1e-7]
    XNear = demand(makeEconomy(shares, endowments, 1 + delta), p);
    verifyLessThan(testCase, max(abs(XNear - XCobbDouglas), [], 'all'), 10 * delta);
    % Utility is continuous in sigma as well.
    uNear = utility(shares(1, :), XCobbDouglas(1, :), 1 + delta);
    verifyEqual(testCase, uNear, utility(shares(1, :), XCobbDouglas(1, :)), 'RelTol', 10 * delta);
end
end

function testFiniteDifferenceStepTradeOff(testCase)
% The central-difference error falls like h^2 and then rises like eps/h,
% so a moderate step beats both a large and a tiny one.
edge = makeEconomy([0.6 0.4; 0.3 0.7], [3 1; 1 3]);
z1 = @(t) [1 0] * excessDemand(edge, [t; 1 - t]);
t = 1.5 / 3.4;
slope = -1.5 / t^2;
centralError = @(h) abs((z1(t + h) - z1(t - h)) / (2 * h) - slope);
verifyLessThan(testCase, centralError(1e-6), centralError(1e-3));
verifyLessThan(testCase, centralError(1e-6), centralError(1e-12));
end

%% Calculus in R^n

function testJacobianIdentities(testCase)
p = [0.7; 1.9; 1.2];
for sigma = [0.5, 1, 2]
    econ = randomEconomy(2, 4, 3, sigma);
    [z, J] = excessDemand(econ, p);
    verifyEqual(testCase, J * p, zeros(3, 1), 'AbsTol', 1e-12);
    verifyEqual(testCase, p' * J, -z', 'AbsTol', 1e-12);
    Jnum = numericalJacobian(@(q) excessDemand(econ, q), p);
    verifyEqual(testCase, J, Jnum, 'AbsTol', 1e-8 * norm(J));
end
end

function testNewtonAgreesWithEigenvector(testCase)
shares = testCase.TestData.shares;
endowments = testCase.TestData.endowments;
cobbDouglas = eigenEquilibrium(makeEconomy(shares, endowments));
p = newtonEquilibrium(makeEconomy(shares, endowments), ones(3, 1));
verifyEqual(testCase, p, cobbDouglas.p, 'RelTol', 1e-12);
pNear = newtonEquilibrium(makeEconomy(shares, endowments, 1 + 1e-6), ones(3, 1));
verifyEqual(testCase, pNear, cobbDouglas.p, 'RelTol', 1e-5);
end

function testNewtonSolvesCesQuadratically(testCase)
shares = testCase.TestData.shares;
endowments = testCase.TestData.endowments;
for sigma = [0.5, 2]
    econ = makeEconomy(shares, endowments, sigma);
    [p, history] = newtonEquilibrium(econ, ones(3, 1));
    verifyEqual(testCase, excessDemand(econ, p), zeros(3, 1), 'AbsTol', 1e-12);
    % Once steps are small, each is at most a constant times the square
    % of the previous one (steps near rounding level are excluded).
    s = history.Step;
    k = find(s(1:end-1) < 0.2 & s(2:end) > 1e-11);
    verifyNotEmpty(testCase, k);
    verifyLessThan(testCase, s(k + 1) ./ s(k).^2, 1);
end
end

function testLagrangeConditionAtDemand(testCase)
shares = testCase.TestData.shares;
endowments = testCase.TestData.endowments;
sigma = 0.5;
p = [1; 2; 0.5];
X = demand(makeEconomy(shares, endowments, sigma), p);
for i = 1:3
    [u, g] = utility(shares(i, :), X(i, :), sigma);
    priceIndex = (shares(i, :) * p.^(1 - sigma))^(1 / (1 - sigma));
    verifyEqual(testCase, X(i, :) * p, endowments(i, :) * p, 'RelTol', 1e-12);
    verifyEqual(testCase, g', p / priceIndex, 'RelTol', 1e-12);
    verifyEqual(testCase, u, endowments(i, :) * p / priceIndex, 'RelTol', 1e-12);
end
end

function testImplicitFunctionTheorem(testCase)
shares = testCase.TestData.shares;
endowments = testCase.TestData.endowments;
sigma = 0.5;
econ = makeEconomy(shares, endowments, sigma);
p = newtonEquilibrium(econ, ones(3, 1));
[~, J] = excessDemand(econ, p);
transfer = zeros(3);
transfer(1, 1) = -1;
transfer(2, 1) = 1;
[X, income] = demand(econ, p);
dzdtau = (X ./ income)' * (transfer * p);
dpIFT = -J(2:end, 2:end) \ dzdtau(2:end);
% Household 2 owns no food, so only tau >= 0 is feasible: use the
% second-order one-sided difference (-3 f(0) + 4 f(h) - f(2h)) / (2h).
h = 1e-5;
pOne = newtonEquilibrium(makeEconomy(shares, endowments + h * transfer, sigma), p);
pTwo = newtonEquilibrium(makeEconomy(shares, endowments + 2 * h * transfer, sigma), p);
dpResolve = (-3 * p(2:end) + 4 * pOne(2:end) - pTwo(2:end)) / (2 * h);
verifyEqual(testCase, dpIFT, dpResolve, 'RelTol', 1e-6);
end

function testTextbookEconomyHasThreeEquilibria(testCase)
r = 2^(8/9) - 2^(1/9);
[z, dz] = textbookExcessDemand([0.5; 1; 2], r);
verifyEqual(testCase, z, zeros(3, 1), 'AbsTol', 1e-14);
verifyEqual(testCase, sum(sign(-dz)), 1);                % index theorem
q = logspace(-3, 3, 20000);
verifyEqual(testCase, nnz(diff(sign(textbookExcessDemand(q, r)))), 3);
% At (r, q) = (7/9, 1) the derivative in q vanishes, and the curve of
% non-trivial equilibria r(q) passes through the same point.
[~, dzCrossing] = textbookExcessDemand(1, 7/9);
verifyEqual(testCase, dzCrossing, 0, 'AbsTol', 1e-14);
qNear = 1 + [-1e-6, 1e-6];
verifyEqual(testCase, (qNear.^(8/9) - qNear.^(1/9)) ./ (qNear - 1), [7/9, 7/9], 'AbsTol', 1e-6);
end

function testScarfInvariantsAndJacobian(testCase)
p0 = [1.2; 1.0; 0.8];
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
[~, P] = ode45(@(t, p) scarfExcessDemand(p), [0 30], p0, opts);
verifyEqual(testCase, sum(P.^2, 2), repmat(sum(p0.^2), size(P, 1), 1), 'AbsTol', 1e-7);
verifyEqual(testCase, prod(P, 2), repmat(prod(p0), size(P, 1), 1), 'AbsTol', 1e-7);
% z is orthogonal to p and to grad(p1 p2 p3), so it is parallel to their cross product.
z = scarfExcessDemand(p0);
verifyEqual(testCase, norm(cross(z, cross(p0, prod(p0) ./ p0))), 0, 'AbsTol', 1e-14);
J = numericalJacobian(@scarfExcessDemand, ones(3, 1));
verifyEqual(testCase, J, -J', 'AbsTol', 1e-9);
verifyEqual(testCase, sort(abs(imag(eig(J)))), [0; sqrt(3)/4; sqrt(3)/4], 'AbsTol', 1e-8);
end

function testScarfSmallCyclePeriod(testCase)
verifyEqual(testCase, scarfCyclePeriod([1.001; 1; 0.999]), 8 * pi / sqrt(3), 'RelTol', 1e-4);
end

function testInvalidInputsAreRejected(testCase)
verifyError(testCase, @() makeEconomy([0.5 0.4], [1 1]), 'makeEconomy:shares');
verifyError(testCase, @() makeEconomy([0.5 0.5], [1 1], -1), 'MATLAB:validators:mustBePositive');
verifyError(testCase, @() eigenEquilibrium(makeEconomy([0.5 0.5], [1 1], 2)), ...
    'eigenEquilibrium:notCobbDouglas');
end

%% Helpers

function econ = randomEconomy(seed, m, n, sigma)
% A reproducible economy with m consumers, n goods and positive data.
rng(seed);
shares = rand(m, n);
shares = shares ./ sum(shares, 2);
econ = makeEconomy(shares, rand(m, n) + 0.1, sigma);
end
