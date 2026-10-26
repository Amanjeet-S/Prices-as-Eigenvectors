# Proofs

The results behind the code, each checked by a test in
[`tests/testEquilibrium.m`](../tests/testEquilibrium.m).

**Notation.** There are $m$ consumers and $n$ goods. Consumer $i$ has weights
$a_{ij} \ge 0$ with $\sum_j a_{ij} = 1$, endowment $e_i$ and income
$m_i = p \cdot e_i$. Supply is $\omega_j = \sum_i e_{ij} > 0$, and excess
demand is $z(p) = \sum_i x_i(p) - \omega$. A price vector $p \gg 0$ is an
equilibrium if $z(p) = 0$.

## 1. Walras' law
$p \cdot z(p) = 0$ for all $p \gg 0$.

*Proof.* Utility is increasing, so every budget binds:
$p \cdot x_i = p \cdot e_i$. Summing over $i$ gives $p \cdot z(p) = 0$. $\blacksquare$

So $z(p) \in p^\perp$, and only $n-1$ market-clearing equations are independent.
*Test:* `testWalrasLawAndHomogeneity`.

## 2. Homogeneity and the Jacobian identities
(a) $z(tp) = z(p)$ for $t > 0$; (b) $J(p)\ p = 0$; (c) $p^{\mathsf T}J(p) = -z(p)^{\mathsf T}$.

*Proof.*
- (a) Scaling all prices leaves every budget set unchanged.
- (b) Differentiate $t \mapsto z(tp)$ at $t = 1$.
- (c) Differentiate $p \cdot z(p) \equiv 0$ with the product rule. $\blacksquare$

At an equilibrium $p^{\ast}$ is a right and a left null vector of $J$, so
$\mathrm{rank}\ J \le n-1$. If the rank equals $n-1$ (a regular equilibrium),
fixing $p_1$ leaves an invertible reduced Jacobian $\tilde J$. To see this, suppose
$\tilde J y = 0$ and let $w = (0, y)$. By (c), row 1 of $J$ is a combination of
the other rows, so $Jw = 0$. Then $w \in \mathrm{span}(p^{\ast})$, and since
$w_1 = 0$, $w = 0$. *Test:* `testJacobianIdentities`.

## 3. Cobb-Douglas equilibrium is an eigenvector
With $v_j = p_j\omega_j$ and $B_{jk} = \sum_i a_{ij}e_{ik}/\omega_k$, $p$ is an
equilibrium if and only if $Bv = v$. Every column of $B$ sums to 1.

*Proof.* Demand is $x_{ij} = a_{ij}m_i/p_j$, so market $j$ clears when
$\sum_i a_{ij}m_i = v_j$. Substituting $m_i = \sum_k (e_{ik}/\omega_k)\ v_k$ gives
$v = Bv$. Also $\sum_j B_{jk} = \sum_i e_{ik}/\omega_k = 1$. $\blacksquare$

So $\mathbf 1^{\mathsf T}B = \mathbf 1^{\mathsf T}$, and 1 is an eigenvalue of $B$,
because $\det(B - I) = \det(B - I)^{\mathsf T}$. Prices are an eigenvector of the
similar matrix $D^{-1}BD$, where $D = \mathrm{diag}(\omega)$.

If every entry of $B$ is positive, Perron-Frobenius (Hubbard & Hubbard §2.7)
gives a positive eigenvector $w$ for a simple dominant eigenvalue $\rho$. Since
$\rho\ \mathbf 1^{\mathsf T}w = \mathbf 1^{\mathsf T}Bw = \mathbf 1^{\mathsf T}w$, we
get $\rho = 1$, so prices are positive and unique up to scale. When the economy
splits into two islands, $B$ is block diagonal and
$\dim\mathrm{null}(B - I) = 2$. Row reduction stays in $\mathbb{Q}$, so the
prices are rational: $p = (1,\ 507/329,\ 174/329)$.

*Tests:* `testEigenvectorClearsMarkets`, `testReducibleEconomyHasNoUniquePrices`.

## 4. Power iteration converges at rate $|\lambda_2|$
If $B = V\Lambda V^{-1}$, $v_{k+1} = Bv_k$ and $\mathbf 1^{\mathsf T}v_0 = 1$, then
$\lVert v_k - v^{\ast}\rVert  \le \mathrm{cond}(V)\ |\lambda_2|^k\ \lVert v_0 - v^{\ast}\rVert$.

*Proof.* Let $e_k = v_k - v^{\ast} = B^k e_0$ and write $e_0 = Vc$. The first row of
$V^{-1}$ is a left eigenvector for the simple eigenvalue $\lambda_1 = 1$, hence a
multiple of $\mathbf 1^{\mathsf T}$. Since $\mathbf 1^{\mathsf T}e_0 = 0$, this gives $c_1 = 0$. Therefore
$\lVert e_k\rVert  = \lVert V\Lambda^k c\rVert  \le \lVert V\rVert \ |\lambda_2|^k\ \lVert V^{-1}\rVert \ \lVert e_0\rVert$. $\blacksquare$

The steps are summable, so $(v_k)$ is Cauchy (Ross §10). In floating point the
error stops at about $10^{-16}$, so the comparison stops at $k = 18$.
*Test:* `testPowerIterationBound`.

## 5. Existence with two goods
With prices $(t, 1-t)$, $S = \sum_i a_{i1}e_{i2} > 0$ and
$R = \sum_i (1-a_{i1})e_{i1} > 0$, the unique equilibrium is $t^{\ast} = S/(S+R)$.

*Proof.* $z_1(t) = S/t - (S+R)$ is continuous on $(0,1)$, tends to $+\infty$ as
$t \to 0^+$ and tends to $-R < 0$ as $t \to 1^-$. The intermediate value theorem
(Ross §18) gives a zero, which is unique because $z_1$ is strictly decreasing. $\blacksquare$

Bisection makes this constructive. The endpoints of the nested intervals are
monotone and bounded, so they converge (Ross §10) to a common limit $c$. By
continuity, $z_1(c) \ge 0$ and $z_1(c) \le 0$, so $z_1(c) = 0$.
*Test:* `testBisectionMatchesClosedForm`.

## 6. The best finite-difference step
Suppose each evaluation of $f$ has rounding error at most $\varepsilon M$, where
$M$ is the size of the terms being combined. By Taylor's theorem the forward
difference has error at most $\tfrac h2|f''| + 2\varepsilon M/h$, and the central
difference at most $\tfrac{h^2}{6}|f'''| + \varepsilon M/h$.

Minimising over $h$ gives $h \sim \sqrt\varepsilon$ for the forward difference
and $h \sim \varepsilon^{1/3}$ for the central one. `numericalJacobian` uses
$h = \varepsilon^{1/3}\max(1, |p_j|)$. *Test:* `testFiniteDifferenceStepTradeOff`.

## 7. Newton converges quadratically
Suppose $f(x^{\ast}) = 0$, $Df(x^{\ast})$ is invertible and $Df$ is Lipschitz with constant
$M$. Taylor's theorem gives $0 = f(x_k) + Df(x_k)(x^{\ast} - x_k) + r_k$ with
$\lVert r_k\rVert  \le \tfrac M2\lVert x_k - x^{\ast}\rVert^2$. Hence $x_{k+1} - x^{\ast} = Df(x_k)^{-1}r_k$, and
the error is squared at each step (see Kantorovich, Hubbard & Hubbard §2.8).

The code iterates in $\theta = \log p$, so by the chain rule
$\partial z_i/\partial\theta_j = J_{ij}\ p_j$. *Test:* `testNewtonSolvesCesQuadratically`.

## 8. CES demand from Lagrange multipliers
Maximise $u(x) = \big(\sum_j a_j^{1/\sigma}x_j^{(\sigma-1)/\sigma}\big)^{\sigma/(\sigma-1)}$
subject to $p \cdot x = m$. The solution is
$$x_j = \frac{a_j p_j^{-\sigma} m}{\sum_k a_k p_k^{1-\sigma}}, \qquad \lambda = \frac1P, \qquad u = \frac mP,$$
where $P = \big(\sum_k a_k p_k^{1-\sigma}\big)^{1/(1-\sigma)}$.

*Proof.* The condition $\nabla u = \lambda p$ gives $x_j = a_j(\lambda p_j)^{-\sigma}u$.
Substituting into $u$ gives $\lambda P = 1$, and the budget constraint gives
$m = uP$. $u$ is homogeneous of degree one and quasi-concave, hence concave, so
this critical point (Hubbard & Hubbard §3.7) is the maximum. $\blacksquare$

Indirect utility is $m/P$, so $\lambda = \partial V/\partial m$ (the envelope
theorem). *Test:* `testLagrangeConditionAtDemand`.

## 9. Comparative statics
Fix $p_1 = 1$ and move endowments by $\tau$ along a transfer, so that
$\sum_i \Delta e_i = 0$. At a regular equilibrium the implicit function theorem
(Hubbard & Hubbard §2.10) gives
$$\tilde p'(0) = -\big[D_{\tilde p}F\big]^{-1}D_\tau F, \qquad
\partial z_j/\partial\tau = \sum_i (x_{ij}/m_i)(\Delta e_i \cdot p).$$

Household 2 owns no food, so only $\tau \ge 0$ is feasible. The check uses the
second-order one-sided difference $(-3\tilde p(0) + 4\tilde p(h) - \tilde p(2h))/(2h)$.
*Test:* `testImplicitFunctionTheorem`.

## 10. Three equilibria and a non-manifold point
In the economy of Mas-Colell, Whinston & Green (§15.B),
$z_2(q, r) = q^{-1/9} - q^{-8/9} + r(1/q - 1)$.

- **The equilibrium set** is the line $q = 1$ together with the curve
  $r = (q^{8/9} - q^{1/9})/(q - 1)$, since $z_2$ is affine in $r$.
- **At $r = 2^{8/9} - 2^{1/9}$** the equilibria are $q = 1/2, 1, 2$.
- **At $(7/9, 1)$** we have $\partial z_2/\partial r = 0$ and
  $\partial z_2/\partial q = 7/9 - r = 0$. The implicit function theorem fails,
  two branches cross, and the set is not a manifold there.
- **Index.** $z_2 \to +\infty$ as $q \to 0$ and $z_2 \to -r$ as $q \to \infty$, so
  the signs of $\partial z_2/\partial q$ at regular zeros alternate and sum
  (negated) to $+1$. The outer equilibria are stable under $\dot q = z_2$; the
  middle one is not.

*Test:* `testTextbookEconomyHasThreeEquilibria`.

## 11. Scarf's cycles
For $z_i = \frac{p_i}{p_i+p_{i+1}} + \frac{p_{i-1}}{p_{i-1}+p_i} - 1$ (indices mod 3),
along $\dot p = z(p)$:

- **Two conserved quantities.** $\frac{d}{dt}\lVert p\rVert^2 = 2p \cdot z = 0$ by Walras' law.
  Also $z_i/p_i = \frac1{p_i+p_{i+1}} - \frac1{p_{i-1}+p_i}$ telescopes, so
  $\frac{d}{dt}\log(p_1p_2p_3) = 0$.
- **Direction of motion.** $z$ is orthogonal to both $p$ and
  $\nabla(p_1p_2p_3)$, so it is parallel to $p \times \nabla(p_1p_2p_3)$.
- **Behaviour near $(1,1,1)$.** There $Jv = \frac14(1,1,1)\times v$, which is
  skew-symmetric, so its eigenvalues are $0$ and $\pm i\sqrt3/4$. Small cycles
  therefore have period $8\pi/\sqrt3$.

Every path stays on the intersection of a sphere with a level set of
$p_1p_2p_3$, a curve that avoids the equilibrium.
*Tests:* `testScarfInvariantsAndJacobian`, `testScarfSmallCyclePeriod`.

## Scope
Existence with many goods needs a fixed-point theorem: Brouwer's when excess
demand is continuous (Mas-Colell, Whinston & Green, Ch. 17), and a
generalisation of Kakutani's in Arrow & Debreu (1954). Both are beyond this
course. Perron-Frobenius, Kantorovich and the implicit function theorem are
cited rather than proved.
