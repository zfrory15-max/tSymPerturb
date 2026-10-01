# Interpretation limits, missing specifications and numerical policy

## Scientific boundaries

The source is a predictive CLPN perturbation framework. Temporal precedence, directed edges, a model operation called an intervention, and optimized scores do not identify causal treatment effects. Stable between-person differences, omitted causes, measurement error and interval choice can affect CLPN coefficients (manuscript pp. 2,11–12). Target anchors require substantive justification and feasibility evidence. Synthetic labels such as Fatigue or Rumination are not clinical target recommendations.

The paper verifies a continuous Gaussian, sparse linear generating model with fixed-penalty lasso. It does not establish general performance for ordinal/zero-inflated data, boundaries, missing-not-at-random dropout, latent confounding, measurement non-invariance, participant heterogeneity, nonlinear/time-varying dynamics, or arbitrary p/n ratios (pp. 11,19–20). Require dedicated validation; do not invent automatic missing-data treatment, ordinal fitting, or person-specific causal inference.

Linear means must be linear in a linked dose, and independent source perturbations must add. Threshold scores, clipping, constraints and nonlinear utility can change these identities only because additional nonlinear operations were introduced. Report the operation explicitly. Neither positive pair incremental value nor positive unsigned communication capacity establishes synergy or signed benefit.

True temporal sequence simulation needs observed multistep transitions or explicit stationarity. Even a stable repeated B does not establish that a real intervention would follow the assumed state update. With actual multi-wave data use wave-specific B unless invariance is supported. Two-wave acquisition ordering is a decision objective only.

## Information absent from the manuscript

The document contains no supplemental appendix, complete numerical B/Sigma1/Psi/mu1 vectors, full module membership list, factor loadings, generation code, exact seed spawning procedure or original source-data tables. A new synthetic demo must be labelled illustrative and MUST NOT claim to reproduce the 22-node population, published top five, .0057-SD Monte Carlo discrepancy, .76/.88/.93 rank recovery, or Fatigue->Pain sequence utility .592.

The manuscript leaves these implementation decisions unspecified:

- Q_H horizon H and propagation discount gamma
- Exact per-target partner sets T_i and their selection procedure
- Numeric outcome weights, clinical utility weights beyond equal-weight methodological verification, and direction-alignment customization
- Sequence dose/update policy u_t, state dependence, constraints, allowable repetitions, cost values, acquisition discount, temporal discount eta and optimizer
- Per-wave standardization ddof, estimator convergence tolerance/iterations, tie ranking and tie-breaking convention
- Degenerate-denominator handling outside the explicit 1e-12 in spillover fraction
- Near-constant min-max tolerance, PSD/symmetry tolerance, floating-point null-test tolerance
- Exact uncertainty interval methods for applied data, bootstrap count, handling of failed fits, and tied top-K membership

Require configuration or label software choices transparently. A default chosen for a runnable example is not a manuscript parameter. Do not silently convert the printed Eq.44 raw-response utility to standardized response; if desired, expose and label a separate convention.

There is a prose nuance on p.9: it reports order-invariant cumulative two-wave acquisition utility without updating/constraints. Equation 45 includes a discount. The mathematically valid invariant is for a fixed final set with no discount (gamma=1), additive U and order-independent costs; discounting can produce acquisition-order differences. Follow the equation and state the assumptions.

## Numerical and input safeguards (engineering policy, not new science)

Reject nonfinite entries, nonsquare/inconsistent B, mismatched node ordering, invalid/duplicate target IDs, invalid module labels, unknown partners, self-pairs, doses or block fractions outside [0,1], noninteger or nonpositive H, and invalid discount intervals. Validate covariance symmetry and positive semidefiniteness within a documented tolerance. Never silently repair a substantially indefinite matrix. Residual covariance may be non-diagonal in the general model; the simulation alone uses independent residuals.

Require positive unperturbed outcome SDs for every standardized outcome. Zero source SD prevents standardization/anchor transformation even though a pre-fitted raw-scale query can algebraically exist. Reject or explicitly omit unsupported outcomes with reported scope changes; never replace zero SD with epsilon without disclosure.

Use nonnegative finite outcome weights with positive sum on EACH outcome comparison set. The printed paper does not explicitly constrain outcome w, but this validation is needed for the usual burden/fraction interpretation. Utility weights omega must be nonnegative with positive total. If all comparison weights are zero, cross/pair means are undefined.

Breadth requires p>1; pair utility excluding two outcomes requires p>2; module reach requires at least two nonempty modules. An empty partner set makes Eq.39 undefined. A full seven-dimensional score should reject these cases or mark unavailable, rather than silently setting zero and presenting a manuscript-complete tVPPS. Partial profiles may be reported as such.

If Q_H(B)=0 (e.g. B=0), Eq.35 is 0/0 and undefined. Any zero-return or omitted-dimension convention must be explicitly labeled a software extension. For finite H there is no mathematical requirement rho(B)<1, but large powers can overflow. Detect nonfinite propagation; do not silently clip it. Infinite-horizon interpretations would need additional convergence assumptions, especially rho(gamma*abs(B))<1, not merely rho(B)<1.

Use an explicit documented policy for numerically near-constant dimensions. A small observed range is not scientifically a constant dimension; overaggressive tolerance can erase meaningful rankings. Preserve signed pair increments even though Eq.39 clips their score contribution to zero. Epsilon in Eq.51 is exactly the manuscript's stabilizer; do not add epsilon to other denominators and pretend the resulting quantity is the same equation.

For bounded-normal expectations require L<=U and s>=0. At s=0 return clip(m,L,U). Avoid 0*infinity if supporting infinite bounds; otherwise validate finite bounds. Bounded outcome means need the correct perturbed covariance and a justified Gaussian model; arbitrary marginal clipping does not create an identified causal effect.

## Automated scope versus analytical workflow

A deterministic NumPy implementation of fitted-model moments and operators does not perform CLPN model selection, measurement validation, subject resampling, bootstrap fitting or causal identification. These steps should be documented and delegated to validated analysis code when needed. Full-pipeline uncertainty must repeat fitting, standardization, anchors, scoring and normalization in each participant bootstrap; a model-only API must not describe itself as performing that bootstrap.

Always report source and outcome times, scale and anchor, B orientation, estimator provenance, dose map, candidate/partner sets, modules, thresholds, weights, block convention, H/discounts, signed response profile, raw seven utilities, normalized scores, uncertainty availability and any stationarity assumption. Preserve a distinction between source manuscript equations, transparent engineering conventions and unimplemented research extensions.
