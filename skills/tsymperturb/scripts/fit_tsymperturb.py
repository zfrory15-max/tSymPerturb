"""Optional paired two-wave Lasso fitting and participant-bootstrap adapter.

This is applied-data uncertainty, NOT the manuscript's independent-dataset
Monte Carlo recovery experiment. See references/fitting-and-uncertainty.md.
"""
from __future__ import annotations
import argparse
from collections import Counter
import json
from pathlib import Path
from statistics import NormalDist
import warnings
import numpy as np
from tsymperturb import Model, score, serializable, transform_anchors


def _pairs(t1, t2):
    x, y = np.asarray(t1, dtype=float), np.asarray(t2, dtype=float)
    if x.ndim != 2 or x.shape != y.shape or x.shape[0] < 3 or x.shape[1] < 3:
        raise ValueError('T1 and T2 must be matched n-by-p arrays with n,p >= 3')
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Nonfinite/missing data unsupported: supply prespecified complete matched observations')
    return x, y


def fit_pairs(t1, t2, *, labels, modules, raw_anchors, higher_is_worse,
              alpha=0.03, max_iter=10000, tol=1e-8):
    """Separate wave z scores (ddof=0), per-outcome fixed-alpha Lasso, no intercept.

    Each row must be the same independent participant at both waves. Symptom
    orientation must be performed upstream and explicitly affirmed by caller.
    """
    from sklearn import __version__ as sklearn_version
    from sklearn.linear_model import Lasso
    from sklearn.exceptions import ConvergenceWarning
    x, y = _pairs(t1, t2)
    if higher_is_worse is not True:
        raise ValueError('Explicit higher_is_worse=true required; orient data and anchors upstream')
    if not np.isfinite(alpha) or alpha <= 0:
        raise ValueError('alpha must be finite and positive')
    if isinstance(max_iter, bool) or not isinstance(max_iter, int) or max_iter < 1:
        raise ValueError('max_iter must be a positive integer')
    if not np.isfinite(tol) or tol <= 0:
        raise ValueError('tol must be finite and positive')
    n, p = x.shape
    if len(labels) != p or len(modules) != p:
        raise ValueError('labels and modules must match data columns')
    center1, center2 = x.mean(0), y.mean(0)
    sd1, sd2 = x.std(0, ddof=0), y.std(0, ddof=0)
    if np.any(sd1 <= 0) or np.any(sd2 <= 0):
        raise ValueError('Constant symptom at T1 or T2: standardization undefined')
    z1, z2 = (x-center1)/sd1, (y-center2)/sd2
    B = np.empty((p, p))
    with warnings.catch_warnings():
        warnings.simplefilter('error', ConvergenceWarning)
        for j in range(p):
            fit = Lasso(alpha=alpha, fit_intercept=False, max_iter=max_iter,
                        tol=tol, selection='cyclic').fit(z1, z2[:, j])
            B[j] = fit.coef_
    residual = z2-z1@B.T
    residual -= residual.mean(0)
    # Empirical moments with the same ddof=0 convention as the z-score transform.
    cov1 = z1.T@z1/n
    psi = residual.T@residual/n
    transformed = transform_anchors(raw_anchors, center1, sd1)
    model = Model(source_labels=labels, outcome_labels=labels, B=B,
                  mu1=z1.mean(0), sigma1=cov1, psi=psi,
                  intercept=np.zeros(p), anchors=transformed, modules=modules,
                  orientation='outcome_by_source', higher_is_worse=True,
                  scale='Separate T1/T2 sample z scores, ddof=0')
    diagnostics = dict(n=n, source_center=center1, source_sd=sd1,
                       outcome_center=center2, outcome_sd=sd2,
                       transformed_anchors=transformed, B=B, psi=psi,
                       alpha=alpha, max_iter=max_iter, tol=tol, ddof=0,
                       fit_intercept=False,
                       normalization_sd=model.sd2,
                       numpy_version=np.__version__, sklearn_version=sklearn_version,
                       empirical_predictor_residual_covariance=z1.T@residual/n)
    return model, diagnostics


def _wilson(successes, total, level):
    z = NormalDist().inv_cdf((1+level)/2)
    phat = successes/total
    denom = 1+z*z/total
    center = (phat+z*z/(2*total))/denom
    radius = z*np.sqrt(phat*(1-phat)/total+z*z/(4*total*total))/denom
    return [0.0 if successes == 0 else float(max(0,center-radius)),
            1.0 if successes == total else float(min(1,center+radius))]


def bootstrap_pairs(t1, t2, *, fitting, scoring, n_boot=1000, seed=0,
                    confidence=0.95):
    """Complete-pipeline iid paired-participant bootstrap; no silent retries.

    Competition ranks: 1 + number of strictly larger scores. Every tied target
    at a top-K boundary is included, so selection probabilities can sum above K.
    """
    x, y = _pairs(t1, t2)
    if isinstance(n_boot, bool) or not isinstance(n_boot, int) or n_boot < 1:
        raise ValueError('n_boot must be a positive integer')
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError('seed must be a nonnegative integer')
    if not np.isfinite(confidence) or not 0 < confidence < 1:
        raise ValueError('confidence must lie in (0,1)')
    original_model, original_fit = fit_pairs(x, y, **fitting)
    point = score(original_model, **scoring)  # Validate full point pipeline first.
    rng = np.random.default_rng(seed)
    replicates, failures = [], []
    for b in range(n_boot):
        indices = rng.integers(0, len(x), size=len(x))
        try:
            model, diagnostics = fit_pairs(x[indices], y[indices], **fitting)
            result = score(model, **scoring)
            replicates.append(dict(attempt=b+1, fit=diagnostics, analysis=result))
        except (ValueError, ArithmeticError, np.linalg.LinAlgError, Warning) as exc:
            failures.append(dict(attempt=b+1, error_type=type(exc).__name__, reason=str(exc)))
    succeeded = len(replicates)
    summary = []
    if succeeded:
        scores = np.asarray([r['analysis']['tvpps'] for r in replicates])
        ranks = np.asarray([r['analysis']['ranks'] for r in replicates])
        q = [(1-confidence)/2, (1+confidence)/2]
        for i, label in enumerate(point['candidates']):
            selection = {}
            for k in (1,3,5):
                count = int(np.sum(ranks[:,i] <= k))
                selection[str(k)] = dict(count=count, probability=count/succeeded,
                                         monte_carlo_wilson_interval=_wilson(count,succeeded,confidence))
            summary.append(dict(target=label,
                                tvpps_percentile_interval=np.quantile(scores[:,i],q),
                                rank_percentile_interval=np.quantile(ranks[:,i],q),
                                top_k=selection))
    return dict(point_estimate=point, point_fit=original_fit, attempted=n_boot,
                succeeded=succeeded, failed=len(failures), seed=seed,
                confidence=confidence, fitting=fitting, scoring=scoring,
                status='no_successful_replicates' if not succeeded else
                       ('partial_failures' if failures else 'ok'),
                failure_reason_counts=dict(Counter(f"{f['error_type']}: {f['reason']}" for f in failures)),
                failures=failures, summary=summary, replicates=replicates,
                tie_rule='Competition rank = 1 + count(strictly larger tVPPS); include all ties at K',
                denominator='Successful replicates only; failures may bias conditional summaries',
                interpretation='Model-implied uncertainty under iid participants; not causal effects or independent-dataset Monte Carlo recovery')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, help='JSON with t1,t2,fitting,scoring')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--n-boot', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=0)
    args = parser.parse_args()
    try:
        data = json.loads(args.input.read_text(encoding='utf-8'))
        result = bootstrap_pairs(**data, n_boot=args.n_boot, seed=args.seed)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(serializable(result), ensure_ascii=False,
                                         allow_nan=False, indent=2)+'\n', encoding='utf-8')
    except (ValueError, TypeError, KeyError, OSError, ImportError, Warning) as exc:
        parser.exit(2, f'Fitting/bootstrap error: {exc}\n')

if __name__ == '__main__':
    main()
