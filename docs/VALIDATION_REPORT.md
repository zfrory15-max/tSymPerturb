# Independent validation report

Initial audit: 2026-09-30. Updated and independently rerun: 2026-10-01. Scope: equation-level implementation and user-facing workflow. This audit is separate from code authoring: expected values were independently derived from the supplied manuscript before evaluating the implementation.

## Status and meaning

**Core numerical audit: PASS.** No unresolved equation-level defect was found in the tested cases. This is not proof for every input, clinical validation, causal identification, or replication of the manuscript's original simulation results.

| Check | Status | Evidence |
|---|---|---|
| All seven raw utilities and weighted tVPPS | PASS | Independent 22-symptom/four-module calculation plus an internal asymmetric 4-node hand-check fixture with unequal weights, nonzero anchors and correlated sources |
| Candidate-set normalization | PASS | Direct equation oracle; exact constant dimension and singleton candidate set produce 50 |
| Matrix direction and signed effects | PASS | Columns are outgoing; edge block changes B[outcome,source]; negative responses retained |
| State and covariance propagation | PASS | Dose endpoints and partial doses; cross-covariance scaling; location-only versus linked map; baseline outcome SD fixed |
| Standardized anchor transformation | PASS | Raw anchors transformed by source mean and source SD |
| Linear null constraints | PASS | Dose grid 0, .001, .2, .7, 1; joint response equals sum of singles; pair calculations exclude both targets |
| Communication blocking | PASS | Full column, cross-only diagonal preservation, reference-state response, centered reference zero; finite matrix-power oracle |
| Module and threshold semantics | PASS | Additional unequal-size-module fixture; unweighted per-module means; inclusive thresholds |
| Multi-wave recursion and acquisition | PASS | Noncommuting matrices; intervention before transition; discounted order and undiscounted fixed-set ties |
| Clipped-normal mean | PASS | Independent numerical quadrature, including zero variance and outside-bound means |
| Invalid inputs | PASS | Orientation, label order, covariance, dose, partner, zero comparison weight, zero propagation capacity, discount and horizon rejected |
| Derived overflow | PASS after repair | Large finite B formerly produced infinite baseline SD; now rejected. Overflowing weight totals also rejected |
| Independent Monte Carlo | PASS | 22 targets × 120,000 common-noise draws; seed 20261002; maximum absolute mean error / Monte Carlo SE = 2.728486; maximum covariance error = 0.004429 raw covariance units |
| Isolated JSON CLI | PASS | Input and output in fresh /tmp directory; cwd=/tmp; output checked against independent oracle |
| Maintainer pytest suite | PASS | 55 passed in the independent rerun, including constructed-22 examples and original-input gate tests |
| Skill instruction forward workflow | PASS | Copied skill to a fresh /tmp install directory and ran the exact documented skill-folder command; seven dimensions and dose outputs returned |
| Optional fitting/bootstrap adapter | PASS | 250 paired participants × 22 symptoms; independent fit audit 8/8 replicates; separate CLI check with --n-boot 10 --seed 42: 10/10 successful, 22 target summaries, no failures |
| Bootstrap failure handling | PASS | Real three-participant degeneracy fixture: 30 attempts, 25 successes, 5 failures recorded without replacement |
| Strategy example | PASS | 22-symptom system; multi-transition propagation, repeated additional improvement, bounds and two-target acquisition among eight prespecified candidates |
| Constructed-22 full check | PASS | All seven utilities, 132 target–dose checks and 231 unique pairs; maximum dose response error 2.220446049250313e−16, pair response error 5.551115123125783e−17 |
| Original 22-node experiment reproduction | NOT RUN / insufficient inputs | Manuscript does not supply full generating matrices, algorithm or source data |
| Clinical validity / causal efficacy | NOT ESTABLISHED | Synthetic model checks cannot establish patient benefit or identify treatment effects |

## Constructed example and internal oracle

The public examples use 22 symptoms, four modules (6/6/5/5) and 250 paired participants. All labels and parameters are constructed, with seed 20261001. The system has 53 off-diagonal paths (47 positive, six negative) and spectral radius 0.7675154569729512. It was not tuned to match the paper's spectral radius or ranking. The independently checked top five in this constructed system are S19, S14, S09, S12 and S17.

`audit_synthetic22.py` derives the seven utilities and normalization directly from the supplied arrays, then compares them with the implementation. It checks all 22 targets at six doses and all 231 unique pairs. These results validate the example workflow, not the manuscript's original generating system.

A small four-node fixture remains solely for transparent internal hand checks. The four-node fixture uses an asymmetric transition with positive and negative entries, source means [2,1,3,2.5], anchors [0,-.5,1,0], outcome weights [1,2,3,4], correlated source covariance, and non-identical residual variances. Propagation H=3 and gamma=.7 are explicit audit choices, not claimed manuscript parameters. All targets are candidates and all other nodes are partners. Expected tVPPS is approximately [53.44691970,26.89355622,60.49736087,38.59445372]. Tests compare raw utilities as well as final scores; a plausible ranking alone is insufficient evidence.

Monte Carlo draws verify the equation in a new illustrative system and are not independent clinical datasets. Common random numbers reduce simulation noise. The manuscript's .0057-SD discrepancy, population ranking and finite-sample recovery statistics are not being reproduced.

## Reproduce the independent checks

From the repository root, with NumPy and SciPy available:

```sh
python docs/validation/run_audit.py
python docs/validation/audit_synthetic22.py
python docs/validation/mc_audit.py
python docs/validation/fit_audit.py  # requires optional fitting dependencies
python -m pytest -q
python examples/strategy_demo.py
```

The equation oracle is `docs/validation/oracle.py`. It computes expected utilities independently and does not call production functions. The forward audit imports the public API and invokes the actual JSON CLI from an isolated directory. The public-data audit uses the same constructed 22-symptom inputs as the documented CLI.

## Defects found and repaired

1. Finite but extremely large transition coefficients could overflow derived covariance and pass Model construction with infinite follow-up SD. The builder added finite-derived-moment checks and arithmetic overflow rejection. Independent rerun passed.
2. Individually finite weights could sum to infinity. The builder added nonfinite-total rejection. Independent rerun passed. This conservative policy rejects such inputs rather than rescaling them silently.
3. A final finite-weight hardening change normalizes utility weights before the composite dot product, avoiding intermediate overflow even when the total is finite; the added regression test and independent reruns passed.
4. The integrated suite exposed Wilson probability-interval roundoff at zero successes. Exact endpoint handling was repaired; the full suite and independent adapter checks now pass. A temporary missing fitting-reference link was also resolved before the final suite.
5. Replacing the public four-symptom example exposed nine stale test assumptions about labels and dimensions. These were updated to the 22-symptom inputs; the independent rerun passed all 55 tests.

## Scope boundaries

- The paper omits exact H, gamma, partner choices and complete generating parameters; example choices are transparent software choices.
- Reusing a two-wave B beyond one transition requires stationarity. A two-wave target-acquisition objective is not biological treatment order.
- Absolute propagation capacity is not signed clinical benefit; positive incremental combination value is not causal synergy.
- Rank uncertainty must be produced by a full paired-participant bootstrap, not by changing a fitted matrix alone.
- Passing tests does not establish robustness to arbitrary missingness, ordinal data, latent confounding, measurement error or nonlinear dynamics.

## Optional adapter review

The independent adapter audit reconstructs every sampled row index from the seed, independently standardizes both waves, calls scikit-learn directly for coefficient expectations, transforms raw anchors anew, derives the model-implied normalization variance, and independently reruns min–max normalization and rank selection frequency. It invokes the adapter CLI from a fresh temporary directory. Its 8-replicate run is a software smoke check, not adequate precision for research uncertainty estimation.

The adapter documents that penalized residuals are not necessarily orthogonal to predictors. It uses the manuscript structural covariance convention B Sigma1 B-transpose + Psi, which need not equal the empirical standardized T2 covariance, and reports the empirical predictor–residual covariance for scrutiny. This is an explicit modeling limitation, not evidence that residual independence has been verified. Percentile intervals have no validated coverage claim; Wilson intervals quantify finite-bootstrap Monte Carlo error only. Failure summaries are conditioned on successful replicates and can be biased. Participant identities must already be correctly paired; this adapter checks array alignment, not identity.

Validation environment / 验证环境: Python 3.12.14; NumPy 2.3.5; SciPy 1.17.0; scikit-learn 1.8.0.
