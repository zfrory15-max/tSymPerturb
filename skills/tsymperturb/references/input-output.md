# Input/output contract and Python API

## Portable execution

Python 3.10+ and NumPy 1.24–2.x. The skill is self-contained: its `scripts/tsymperturb.py` can be executed directly or imported by adding that scripts directory to Python's module search path. No editable installation or repository-relative import is required.

Repository-root command:

```bash
python skills/tsymperturb/scripts/tsymperturb.py examples/synthetic_model.json --output outputs/synthetic_results.json
```

Installed skill command (with an absolute input/output path):

```bash
python "$CODEX_HOME/skills/tsymperturb/scripts/tsymperturb.py" INPUT.json --output OUTPUT.json
```

`--output` is required and overwrites that path if it exists. Choose a new path when preservation matters. Inputs are not modified. Invalid inputs exit 2 with a diagnostic; no result should be interpreted after failure.

## JSON fields

Top-level keys: `model`, `scoring`, optional `dose_requests`.

### model: Model keyword arguments

All fields are required; numeric arrays must be finite.

- `source_labels`, `outcome_labels`: same unique, nonempty string labels in exactly the same order. Shape p, p>=3 (the complete score needs pair-excluded outcomes).
- `B`: p×p coefficients, row=outcome, column=source. Diagonal retains autoregression.
- `orientation`: literal `outcome_by_source`.
- `higher_is_worse`: literal true, after any necessary scale reversal outside this script.
- `scale`: nonempty description of the fitted scale and transformations.
- `mu1`: p source means, in the fitted source scale.
- `sigma1`: p×p source covariance; symmetric positive semidefinite.
- `psi`: p×p residual covariance; symmetric positive semidefinite. The model assumes source/residual independence in covariance.
- `intercept`: p fitted intercepts, in outcome scale.
- `anchors`: p explicit source-state anchors in the fitted source scale. Use `transform_anchors(raw_anchors, source_center, source_sd)` when standardized.
- `modules`: p nonempty module strings, at least two distinct modules.

All coefficient, moment, intercept and anchor values must refer to compatible fitted coordinates. For standardized X and Y, do not mix raw X moments into a standardized B. If estimating a model externally, record centering/scaling choices and transform anchors anew within each resample.

Derived baseline values are mu2=a+B mu1 and Sigma2=B Sigma1 Bᵀ+Psi; every derived value must be finite and each outcome SD positive. The script does not silently substitute empirical follow-up SD for this model-implied SD.

### scoring: score keyword arguments

Required:

- `candidates`: unique ordered target labels; normalization only uses this set (one candidate yields all neutral dimensions).
- `partners`: map from every candidate to a nonempty list of unique partner labels; a target cannot partner with itself. Partners can be outside the candidate set but must be model symptoms. This is an explicit scientific choice, not automatically changed when candidates change.
- `outcome_weights`: p nonnegative weights with finite positive sum. Every cross-target/pair-excluded comparison set must also retain positive weight.
- `utility_weights`: seven nonnegative weights with finite positive sum, in the dimension order below.
- `horizon`: positive integer H for unsigned communication capacity Q_H.
- `gamma`: 0<gamma<1 for Q_H.

Optional manuscript reference values:

- `q`: communication block fraction, default .8; 0<=q<=1.
- `breadth_threshold`: default .05 SD, >=0.
- `module_threshold`: default .03 SD, >=0.

H=3 and gamma=.7 in the synthetic example are illustrative choices, NOT recovered manuscript settings. Zero Q_H rejects the full composite because its communication ratio is undefined. A zero utility weight does not bypass definition checks for that dimension: this API always returns a complete seven-dimension profile.

### dose_requests

A list of `Model.state` keyword objects, such as:

```json
[{"doses":{"S1":0.5}}, {"doses":{"S2":1.0},"scale_factors":{"S2":1.0}}]
```

Each dose lies in [0,1]. Default residual multiplier is 1-dose (linked location–scale map). Explicit `scale_factors` must name exactly the same target keys; values in [0,1] specify residual multipliers independently. A value of 1 is location-only, so a dose-one location-only intervention is NOT exact zero-variance t-vKO. Non-target coordinates remain unchanged. The API covers attenuating maps, not arbitrary variance-inflating location–scale maps.

## Output

A JSON object containing:

- `candidates`, `dimensions`: labels in output order.
- `raw`, `normalized`: candidate×7 values. Dimensions in order: downstream efficacy, cross-symptom spillover, breadth, cross-module reach, communication-block value, combination value, spillover fraction.
- `tvpps`: candidate composite values on [0,100]; candidate-set-relative only.
- `ranks`: descending competition rank; exact ties share a rank (1,1,3). Floating-point near ties are not merged; inspect numerical differences and uncertainty.
- `standardized_responses`: candidate×outcome response, standardized by fixed baseline model-implied SD.
- `pairs`: each requested directed target/partner pair; comparison-set single utilities, joint utility, signed increment, additive contrast. Reciprocal pairs may be listed twice when requested in both partner sets.
- `config`: weights, horizon, discount, thresholds, q and partner sets, preserving reproducibility.
- `baseline`: model-implied follow-up mean, covariance, SD.
- `scale`, `source_labels`, `interpretation`.
- `dose_results`: original request plus source/outcome mean/covariance and raw/standardized response.

No automatic confidence interval, fitted coefficient, missing-data repair or clinical recommendation is fabricated. Optional fitting/uncertainty lives in its separately documented adapter.

## Python-only operations

Import `tsymperturb` from the script directory.

- `block(model, source, q, outcome=None, cross_only=False)`: blocked B. Source-only blocks the full outgoing column; cross_only preserves Bii. Supplying outcome blocks one edge.
- `reference_block_response(model, blocked_B, reference_state)`: (B−B*)xref; specify reference state in fitted coordinates. This is a coefficient-blocking estimand with the intercept held fixed, not an anchor-subtracted state intervention.
- `propagation_capacity(B,horizon,gamma)`: sum of all entries of sum((gamma abs(B))^h), h=1..H; finite horizons do not require spectral radius<1.
- `bounded_normal_mean(mean,sd,lower,upper)`: Eq28 marginal clipped Gaussian expectation. SD=0 returns clipped mean. It does not derive a joint bounded covariance.
- `trajectories(transitions,initial_improvement,interventions=None)`: returns H×p raw differences; ordered explicit transition list. Recursion delta_next=B_t(delta_current+u_t); `interventions` is H×p additional pre-transition improvement, not a dose or fresh full knockout each time. Zero interventions give Eq41. To model anchoring on a changing state, compute valid u_t yourself from that state and document the rule.
- `horizon_utility(responses,outcome_weights,eta)`: Eq44 on the supplied H×p scale; 0<eta<=1. No implicit SD division.
- `acquisition_order(single_utilities,costs,length,gamma)`: exact exhaustive fixed-length additive Eq45 objective; dictionary keys are target names, costs nonnegative, 0<gamma<=1, at most nine targets. Returns tied best orders (absolute tie tolerance 1e-12). No budget/feasibility or nonlinear-set utility solver is implemented. This is a decision ordering, never identified treatment timing.

Stable infinite-horizon dynamics, nonlinear transitions, state-dependent dosing, changing measurement sets and generalized sequence optimization are outside these helpers. See limitations before extending them.
