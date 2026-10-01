---
name: tsymperturb
description: Apply, audit, or explain tSymPerturb temporal virtual perturbations in longitudinal cross-lagged symptom networks, including t-vKO, t-vKD, dosage, directed blocking, seven-outcome tVPPS, combinations and multi-wave strategies. Use for tSymPerturb or time-indexed CLPN intervention hypotheses, not cross-sectional SymPerturb or causal treatment recommendations.
---

# tSymPerturb

Turn a fitted, directed longitudinal symptom model into auditable model-implied intervention hypotheses. Use the supplied numerical implementation rather than recreating equations from memory.

## Establish the estimand

Read [methods](references/methods.md) for the requested operators and [input/output](references/input-output.md) before calculating. Check source/outcome occasion, measurement interval, symptom order, scale, anchors, candidate and partner sets, modules, burden weights and utility weights. B[j,i] means source i → outcome j; do not transpose by guesswork. Higher values must mean worse states.

Use original-scale clinically meaningful anchors transformed into the fitted scale. A standardized zero is the mean, generally not symptom absence. Preserve the unperturbed outcome SD for all dose comparisons. Never refit a network after hypothetical knockout.

## Calculate

- Run the portable JSON analysis with `python scripts/tsymperturb.py INPUT.json --output OUTPUT.json` from this skill folder. Install NumPy in the execution environment if unavailable. The input reference documents every required field.
- Use `Model.state` for linked t-vKO/t-vKD, or explicitly specified separate scale attenuation; retain predicted covariance as well as mean response.
- Use `block` and `reference_block_response` for signed edge/source-node effects at an explicit reference state. `score` computes unsigned finite-horizon communication capacity separately from signed benefit.
- Use `score` for all seven raw dimensions, candidate-set normalization and tVPPS. Require explicit H, gamma, partner sets and weights. The manuscript supplies q=.8, breadth=.05 SD and module threshold=.03 SD; H/gamma have no published numeric default in the supplied version.
- Use `bounded_normal_mean` only when Gaussian marginal clipping is a declared extra operation; label any curvature as boundary-induced.
- Use `trajectories` with explicit ordered transition matrices and additional pre-transition improvement vectors. Repeating one fitted B beyond its observed interval requires an explicit stationarity assumption. `horizon_utility` uses supplied response units; equation 44 uses raw R.
- `acquisition_order` optimizes only the additive fixed-length two-wave target-acquisition objective (up to nine candidates). It is not a biological treatment-order optimizer. State-dependent temporal interventions need a separately specified update rule.

## Verify and report

Read [limitations](references/limitations.md) before interpretation. Check exact dose linearity and independent-source additivity on a common outcome set. Pair incremental values must exclude both targets, retain signed increments and label the positive-part average as combination value, never synergy.

Report seven raw dimensions beside normalized scores, candidate-set definition, all parameter choices, source/outcome scales and the predicted adverse as well as beneficial responses. Do not treat close ranks as certain. A constant dimension scores 50; undefined denominators are errors, not zeros. Clinical utility weights must be prespecified; equal weights are only a methodological choice.

The base CLI starts from fitted moments. For optional continuous complete-case lasso fitting and participant-level complete-pipeline bootstrap, read [fitting and uncertainty](references/fitting-and-uncertainty.md) and use its separate adapter. Do not call a fixed-network perturbation-only bootstrap full-pipeline uncertainty. Report failures and effective successful replicates.

Distinguish internal equation verification from external/clinical validation. The included synthetic example is illustrative and does not recreate the paper's original 22-node simulation. Never invent a DOI, author authorization, causality, validated treatment policy, or empirical information about repeated treatment order from two waves.
