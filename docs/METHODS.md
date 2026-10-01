# Method overview

[English](METHODS.md) | [简体中文](zh-CN/METHODS.md)

This guide summarizes the supplied manuscript, *tSymPerturb converts longitudinal symptom networks into time-indexed intervention strategies*. The implementation is a reference translation of that specification. It does not establish a causal treatment effect or reproduce the manuscript's original simulation tables.

## 1. Time, direction and scale

For the two-wave model `X2 = a + B X1 + ε`, rows of `B` are T2 outcomes and columns are T1 sources. Thus `B[j,i]` is `i(T1) → j(T2)`; diagonal entries represent autoregressive persistence. All symptoms must be oriented so that larger values mean worse states.

Define anchors on the original measurement scale. When a source is standardized as `(X1[i] − mean1[i]) / sd1[i]`, transform its anchor with the same mean and SD. A standardized zero usually means the sample mean, not symptom absence. Outcome improvements use the prespecified, unperturbed T2 SD, not a perturbed SD. The provided core computes this denominator from the model-implied baseline covariance `B Σ1 Bᵀ + Ψ`; an empirical-SD convention would be a separately documented choice.

## 2. Source-state operators

A temporal virtual knockdown at dose `d` uses the linked map `X1*[i] = c[i] + (1−d)(X1[i]−c[i])`. A temporal virtual knockout is its `d=1` endpoint. Non-target source coordinates stay unchanged. Target variance is multiplied by `(1−d)^2` and target–non-target covariance by `(1−d)`.

The mean response is `R = B(mean1 − mean1*)`. For one target, `R_i(d) = d B[:,i](mean1[i]−c[i])`. Keep the fitted transition matrix fixed during this query; refitting after each hypothetical intervention changes the estimand. Exact knockout creates a constant source coordinate but does not require standardizing that coordinate again.

The covariance response, when required, is `Σ2* = B Σ1* Bᵀ + Ψ`. Mean-only and linked location–scale interventions can share a mean response while producing different follow-up variances.

## 3. Transition operators

Directed edge blocking attenuates one `B[j,i]`; source-node blocking attenuates column `i`. Specify whether the autoregressive path is included. The signed change in prediction is `(B−B*) x_ref`, with intercept held fixed. Its magnitude depends on the reference profile. A mean-centered source population has zero mean signed effect under this calculation.

The manuscript separately defines unsigned route capacity:

`Q_H(B) = 1ᵀ [Σ(h=1…H) (γ |B|)^h] 1`

The relative loss after blocking is a topology diagnostic. It cannot be read as signed symptom benefit because negative and positive paths enter through absolute values. State the horizon, damping and block fraction explicitly; these are analysis choices.

## 4. Seven tVPPS dimensions

| Dimension | Interpretation |
|---|---|
| Downstream efficacy | Weighted standardized benefit across all T2 symptoms, including the target |
| Cross-symptom spillover | Weighted benefit excluding the target's own T2 outcome |
| Breadth | Fraction of non-target outcomes meeting a prespecified benefit threshold |
| Cross-module reach | Fraction of other modules whose mean benefit meets a prespecified threshold |
| Communication-block value | Relative reduction in unsigned finite-horizon route capacity |
| Combination value | Average positive incremental pair value over a prespecified partner set |
| Spillover fraction | Positive non-target benefit divided by total positive benefit |

Each dimension is direction-aligned and min–max normalized within the prespecified candidate set to 0–100. A constant dimension receives 50. tVPPS is a non-negatively weighted average with at least one positive weight. Equal weights are a methodological default, not validated clinical preference weights. Report raw utilities alongside normalized scores.

The temporal dimensions differ from those of cross-sectional SymPerturb. Dose efficiency and low-dose responsiveness are redundant with efficacy in this linear reference model. Robustness and sampling uncertainty remain separate from tVPPS.

## 5. Three falsification checks

1. **Linear dose response:** `G_i(d) = d G_i(1)` when transition, dose map and mean utility are linear and unbounded.
2. **Additivity:** independent source perturbations satisfy `R_{i,k} = R_i + R_k`; a common-outcome linear utility has additive interaction contrast zero.
3. **Temporal identification:** one observed T1→T2 transition does not identify repeated treatment order.

For combination value, use the same outcome set excluding both pair members for the pair and both singles. Retain the signed increment `I_ik = G_ik − max(G_i,G_k)` even though the reference target-level utility averages `max(0,I_ik)`. Incremental value is distinct from synergy.

Bounded outcomes, nonlinear state maps, interactions or state-dependent transitions may break the first two identities; any such extension needs explicit specification and separate validation.

## 6. Multiple waves and sequences

For a one-time perturbation, propagate through observed, ordered transitions: `R_{t+h} = B_{t+h−1} … B_t δ_t`. Reusing one matrix as `B^h δ` beyond the observed interval introduces a stationarity assumption.

For repeated interventions, use `Δ⁻_{t+1} = B_t(Δ⁻_t + u_t)`. Describe intervention timing, horizon, costs, discounting and feasibility. A two-wave target-acquisition ordering is a decision procedure; it is not evidence about biological treatment order.

## 7. Uncertainty and reporting

Applied uncertainty should rerun participant resampling, standardization, anchor transformation, network estimation, perturbation, utility calculation, normalization and ranking. Report rank distributions and top-K selection probabilities separately. Resampling only the final score table is not a complete-pipeline bootstrap. The [optional adapter](../skills/tsymperturb/references/fitting-and-uncertainty.md) implements this pipeline for independent complete matched two-wave observations with fixed-alpha Lasso. It records failures and source–residual covariance; penalized residuals need not be empirically orthogonal, so model-implied and empirical follow-up variances may differ.

At minimum, report waves and intervals; sample and missingness handling; symptom direction, scales and anchors; matrix orientation and estimator; outcome SDs; candidate/partner sets; modules, weights and thresholds; blocking and horizon settings; raw utilities and tVPPS; algebraic checks; and limitations. Distinguish what the current scripts compute from analyses that require additional code and validation.
