"""Auditable linear tSymPerturb operators. B[outcome, source]; no model refitting."""
from __future__ import annotations
import argparse
import itertools
import json
import math
from pathlib import Path
import numpy as np

DIMENSIONS = ('downstream_efficacy', 'cross_symptom_spillover', 'breadth',
              'cross_module_reach', 'communication_block_value',
              'combination_value', 'spillover_fraction')


def array(x, shape=None, name='array'):
    a = np.asarray(x, dtype=float)
    if (shape is not None and a.shape != shape) or not np.all(np.isfinite(a)):
        raise ValueError(f'{name} must be finite with shape {shape}; got {a.shape}')
    return a


def covariance(x, p, name):
    a = array(x, (p, p), name)
    tol = 1e-10 * max(1., np.linalg.norm(a, 2))
    if not np.allclose(a, a.T, atol=tol, rtol=0) or np.linalg.eigvalsh(a).min() < -tol:
        raise ValueError(f'{name} must be symmetric positive semidefinite')
    return (a + a.T) / 2


def fraction(x, name='fraction', strict=False):
    x = float(x)
    if not np.isfinite(x) or not (0 < x < 1 if strict else 0 <= x <= 1):
        raise ValueError(f'{name} must be in {"(0,1)" if strict else "[0,1]"}')
    return x


def weights(x, p, name='weights'):
    a = array(x, (p,), name)
    with np.errstate(over='ignore', invalid='ignore'):
        total = a.sum()
    if np.any(a < 0) or not np.isfinite(total) or total <= 0:
        raise ValueError(f'{name} must be nonnegative and have finite positive sum')
    return a


def labels(x, name):
    if not isinstance(x, (list, tuple)) or not x or any(not isinstance(v, str) or not v.strip() for v in x) or len(set(x)) != len(x):
        raise ValueError(f'{name} must be nonempty unique string labels')
    return tuple(x)


class Model:
    """Moments and fitted directed transition, all in one explicitly declared scale."""
    def __init__(self, *, source_labels, outcome_labels, B, mu1, sigma1, psi,
                 intercept, anchors, modules, orientation, higher_is_worse, scale):
        self.labels = labels(source_labels, 'source_labels')
        if tuple(outcome_labels) != self.labels:
            raise ValueError('source/outcome labels must have identical order; reorder explicitly')
        if orientation != 'outcome_by_source' or higher_is_worse is not True:
            raise ValueError('require outcome_by_source orientation and higher_is_worse=true')
        if not isinstance(scale, str) or not scale.strip():
            raise ValueError('explicit scale description required')
        self.scale = scale
        self.p = p = len(self.labels)
        if p < 3:
            raise ValueError('at least three symptoms required for pair-excluded utilities')
        self.B = array(B, (p, p), 'B')
        self.mu1 = array(mu1, (p,), 'mu1')
        self.sigma1 = covariance(sigma1, p, 'sigma1')
        self.psi = covariance(psi, p, 'psi')
        self.intercept = array(intercept, (p,), 'intercept')
        self.anchors = array(anchors, (p,), 'anchors')
        if len(modules) != p or any(not isinstance(m, str) or not m.strip() for m in modules) or len(set(modules)) < 2:
            raise ValueError('modules must assign every symptom to at least two named modules')
        self.modules = tuple(modules)
        with np.errstate(over='raise', invalid='raise'):
            self.mu2 = array(self.intercept + self.B @ self.mu1, (p,), 'derived mu2')
            self.sigma2 = array(self.B @ self.sigma1 @ self.B.T + self.psi, (p,p), 'derived sigma2')
            self.sd2 = array(np.sqrt(np.diag(self.sigma2)), (p,), 'derived SD2')
        if np.any(self.sd2 <= 0):
            raise ValueError('unperturbed follow-up SD must be positive for every outcome')

    def index(self, label):
        if label not in self.labels:
            raise ValueError(f'unknown symptom: {label}')
        return self.labels.index(label)

    def state(self, doses, scale_factors=None):
        """Eq13–25. doses map labels to mean reduction in [0,1].

        Default linked scale map is 1-dose. Explicit scale_factors map precisely
        the same targets to residual multipliers in [0,1] (1 = location only).
        """
        if not isinstance(doses, dict) or not doses:
            raise ValueError('doses must be a nonempty label-to-dose mapping')
        if scale_factors is not None and set(scale_factors) != set(doses):
            raise ValueError('scale_factors must name exactly the dose targets')
        mu = self.mu1.copy()
        f = np.ones(self.p)
        for label, dose in doses.items():
            i = self.index(label)
            d = fraction(dose, 'dose')
            mu[i] = self.anchors[i] + (1-d) * (self.mu1[i]-self.anchors[i])
            f[i] = 1-d if scale_factors is None else fraction(scale_factors[label], 'scale factor')
        with np.errstate(over='raise', invalid='raise', divide='raise'):
            sig = array(f[:, None] * self.sigma1 * f[None, :], (self.p,self.p), 'perturbed source covariance')
            response = array(self.B @ (self.mu1-mu), (self.p,), 'response')
            out_mean = array(self.intercept+self.B@mu, (self.p,), 'perturbed outcome mean')
            out_cov = array(self.B@sig@self.B.T+self.psi, (self.p,self.p), 'perturbed outcome covariance')
            standardized = array(response/self.sd2, (self.p,), 'standardized response')
        return dict(source_mean=mu, source_covariance=sig, outcome_mean=out_mean,
                    outcome_covariance=out_cov, response=response, standardized_response=standardized)


def transform_anchors(raw_anchors, source_center, source_sd):
    c = array(raw_anchors, name='raw anchors')
    center = array(source_center, c.shape, 'source center')
    sd = array(source_sd, c.shape, 'source SD')
    if c.ndim != 1 or np.any(sd <= 0):
        raise ValueError('anchor vectors must be one-dimensional and SD positive')
    return (c-center)/sd


def block(model, source, q, outcome=None, cross_only=False):
    """Eq29/32/33; full outgoing column by default, never incoming row."""
    q = fraction(q, 'block fraction')
    i = model.index(source)
    b = model.B.copy()
    if outcome is not None:
        if cross_only:
            raise ValueError('edge block cannot also request cross_only')
        b[model.index(outcome), i] *= 1-q
    else:
        b[:, i] *= 1-q
        if cross_only:
            b[i, i] = model.B[i, i]
    return b


def reference_block_response(model, blocked_B, reference_state):
    return (model.B-array(blocked_B, (model.p, model.p), 'blocked B')) @ array(reference_state, (model.p,), 'reference state')


def propagation_capacity(B, horizon, gamma):
    b = array(B, name='B')
    if b.ndim != 2 or b.shape[0] != b.shape[1]:
        raise ValueError('B must be square')
    if isinstance(horizon, bool) or not isinstance(horizon, int) or horizon < 1:
        raise ValueError('horizon must be a positive integer')
    gamma = fraction(gamma, 'gamma', strict=True)
    a = gamma*np.abs(b)
    power = np.eye(len(b))
    total = 0.
    with np.errstate(over='raise', invalid='raise'):
        for _ in range(horizon):
            power = power @ a
            total += power.sum()
    if not np.isfinite(total):
        raise ValueError('propagation capacity overflow')
    return float(total)


def weighted_mean(values, w, mask=None):
    if mask is not None:
        values, w = values[mask], w[mask]
    if w.sum() <= 0:
        raise ValueError('comparison outcome set must have positive weight')
    with np.errstate(over='raise', invalid='raise', divide='raise'):
        result = float(np.dot(values, w/w.sum()))
    if not np.isfinite(result):
        raise ValueError('weighted mean overflow')
    return result


def score(model, *, candidates, partners, outcome_weights, utility_weights,
          horizon, gamma, q=0.8, breadth_threshold=0.05, module_threshold=0.03):
    """Eq35,38–39,46–53. All choices recorded; partners explicitly prespecified."""
    candidates = labels(candidates, 'candidates')
    ids = [model.index(c) for c in candidates]
    if set(partners) != set(candidates):
        raise ValueError('partners must contain exactly every candidate key')
    w = weights(outcome_weights, model.p, 'outcome_weights')
    uw = weights(utility_weights, 7, 'utility_weights')
    for t in (breadth_threshold, module_threshold):
        if not np.isfinite(t) or t < 0:
            raise ValueError('thresholds must be finite and nonnegative')
    base_capacity = propagation_capacity(model.B, horizon, gamma)
    if base_capacity <= 0:
        raise ValueError('Q_H(B)=0: communication ratio undefined; do not fabricate tVPPS')
    with np.errstate(over='raise', invalid='raise', divide='raise'):
        response = array(model.B * (model.mu1-model.anchors)[None, :] / model.sd2[:, None],
                         (model.p,model.p), 'standardized response matrix')
    raw, pairs = [], []
    for name, i in zip(candidates, ids):
        r = response[:, i]
        nonself = np.arange(model.p) != i
        othermods = sorted(set(model.modules)-{model.modules[i]})
        reach = np.mean([r[np.array(model.modules) == m].mean() >= module_threshold for m in othermods])
        partner_names = labels(partners[name], f'partners[{name}]')
        incremental = []
        for other in partner_names:
            k = model.index(other)
            if k == i:
                raise ValueError('target cannot partner with itself')
            mask = (np.arange(model.p) != i) & (np.arange(model.p) != k)
            gi = weighted_mean(r, w, mask)
            gk = weighted_mean(response[:, k], w, mask)
            joint = weighted_mean(r+response[:, k], w, mask)
            inc = joint-max(gi, gk)
            incremental.append(max(0., inc))
            pairs.append(dict(target=name, partner=other, single_target=gi, single_partner=gk,
                              joint=joint, signed_increment=inc, additive_contrast=joint-gi-gk))
        comm = (base_capacity-propagation_capacity(block(model, name, q), horizon, gamma))/base_capacity
        positive = np.maximum(r, 0)*w
        raw.append([weighted_mean(r,w), weighted_mean(r,w,nonself),
                    float(np.mean(r[nonself] >= breadth_threshold)), float(reach), comm,
                    float(np.mean(incremental)), float(positive[nonself].sum()/(positive.sum()+1e-12))])
    raw = array(raw, (len(candidates),7), 'raw utilities')
    lo, hi = raw.min(axis=0), raw.max(axis=0)
    spread = hi-lo
    normalized = np.full_like(raw, 50.)
    varying = spread != 0  # Exact constant dimensions, as specified in Eq52.
    normalized[:, varying] = 100*(raw[:, varying]-lo[varying])/spread[varying]
    composite = array(normalized @ (uw/uw.sum()), (len(candidates),), 'tVPPS')
    # Competition ranks retain exact ties; alphabetical order is display only.
    ranks = [1+int(np.sum(composite > v)) for v in composite]
    return dict(candidates=list(candidates), dimensions=list(DIMENSIONS), raw=raw,
                normalized=normalized, tvpps=composite, ranks=ranks,
                standardized_responses=response[:, ids].T, pairs=pairs,
                config=dict(outcome_weights=w, utility_weights=uw, horizon=horizon, gamma=gamma,
                            q=q, breadth_threshold=breadth_threshold, module_threshold=module_threshold,
                            partners=partners),
                interpretation='Model-implied, candidate-set-relative hypotheses; not causal treatment effects.')


def bounded_normal_mean(mean, sd, lower, upper):
    """Eq28 for independent marginal clipping, including sd=0; not re-fitting."""
    m,s,l,u = np.broadcast_arrays(*[array(x) for x in (mean,sd,lower,upper)])
    if np.any(s < 0) or np.any(l >= u):
        raise ValueError('require SD >=0 and lower < upper')
    safe = np.where(s == 0, 1, s)
    a,b = (l-m)/safe, (u-m)/safe
    cdf = np.vectorize(lambda x: .5*(1+math.erf(x/math.sqrt(2))), otypes=[float])
    phi = lambda x: np.exp(-x*x/2)/math.sqrt(2*math.pi)
    result = l*cdf(a)+m*(cdf(b)-cdf(a))+s*(phi(a)-phi(b))+u*(1-cdf(b))
    return np.where(s == 0, np.clip(m,l,u), result)


def trajectories(transitions, initial_improvement, interventions=None):
    """Eq41/43; explicit transition list; u is additional pre-transition improvement."""
    delta = array(initial_improvement, name='initial improvement')
    if delta.ndim != 1 or not transitions:
        raise ValueError('require vector improvement and nonempty transition list')
    p = len(delta)
    bs = [array(b, (p,p), 'transition') for b in transitions]
    us = np.zeros((len(bs),p)) if interventions is None else array(interventions, (len(bs),p), 'interventions')
    result = []
    for b,u in zip(bs,us):
        delta = b @ (delta+u)
        if not np.all(np.isfinite(delta)):
            raise ValueError('trajectory overflow')
        result.append(delta.copy())
    return np.asarray(result)


def horizon_utility(responses, outcome_weights, eta):
    """Eq44 in supplied response units (raw R in manuscript), no hidden SD scaling."""
    r = array(responses, name='responses')
    if r.ndim != 2 or not len(r):
        raise ValueError('responses must be nonempty horizon by outcome matrix')
    w = weights(outcome_weights, r.shape[1])
    eta = fraction(eta, 'eta')
    if eta == 0:
        raise ValueError('eta must be >0')
    with np.errstate(over='raise', invalid='raise'):
        result = float(np.dot(eta**np.arange(len(r)), r@(w/w.sum())))
    if not np.isfinite(result):
        raise ValueError('horizon utility overflow')
    return result


def acquisition_order(single_utilities, costs, length, gamma):
    """Exact Eq45 for ADDITIVE set utility only; not biological temporal order.

    Searches all fixed-length unique orders, max nine candidates to bound factorial work.
    Returns all tied best orders (within 1e-12) and their common objective.
    """
    names = labels(list(single_utilities), 'targets')
    if set(costs) != set(names) or len(names) > 9:
        raise ValueError('costs must match targets; exact search supports at most nine')
    u = array([single_utilities[n] for n in names], (len(names),), 'utilities')
    c = array([costs[n] for n in names], (len(names),), 'costs')
    if np.any(c < 0) or isinstance(length,bool) or not isinstance(length,int) or not 1 <= length <= len(names):
        raise ValueError('nonnegative costs and valid positive length required')
    gamma = fraction(gamma, 'acquisition gamma')
    if gamma == 0:
        raise ValueError('acquisition gamma must be >0')
    best, orders = -np.inf, []
    for order in itertools.permutations(range(len(names)), length):
        val = sum(gamma**m*u[i]-c[i] for m,i in enumerate(order))
        if val > best+1e-12:
            best,orders = val,[order]
        elif abs(val-best) <= 1e-12:
            orders.append(order)
    return dict(objective=float(best), orders=[[names[i] for i in o] for o in orders],
                interpretation='Additive target-acquisition decision only; not treatment timing.')


def serializable(x):
    if isinstance(x, np.ndarray): return x.tolist()
    if isinstance(x, np.generic): return x.item()
    if isinstance(x, dict): return {k:serializable(v) for k,v in x.items()}
    if isinstance(x, (tuple,list)): return [serializable(v) for v in x]
    return x


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, help='JSON with model, scoring and optional dose_requests')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        data = json.loads(args.input.read_text(encoding='utf-8'))
        model = Model(**data['model'])
        result = score(model, **data['scoring'])
        result['baseline'] = dict(mean=model.mu2, covariance=model.sigma2, sd=model.sd2)
        result['scale'] = model.scale
        result['source_labels'] = list(model.labels)
        result['dose_results'] = [dict(request=r, result=model.state(**r)) for r in data.get('dose_requests', [])]
        text = json.dumps(serializable(result), ensure_ascii=False, indent=2, allow_nan=False)+'\n'
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding='utf-8')
    except (ValueError, TypeError, KeyError, OSError, FloatingPointError) as exc:
        parser.exit(2, f'Input/analysis error: {exc}\n')

if __name__ == '__main__':
    main()
