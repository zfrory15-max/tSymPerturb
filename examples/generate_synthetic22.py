"""Generate the NEW 22-symptom illustrative example, never manuscript recovery.

Reproducibility: NumPy PCG64 SeedSequence(20261001), separate graph/sample streams.
All names, module assignments, edge choices, factor loadings, scoring settings,
and samples are constructed here. Matching a few disclosed design constraints
cannot authenticate this system as the unavailable original generating system.
"""
import argparse
import json
from pathlib import Path
import numpy as np

SEED = 20261001
P = 22
N = 250


def generate():
    graph_seed, sample_seed = np.random.SeedSequence(SEED).spawn(2)
    graph = np.random.Generator(np.random.PCG64(graph_seed))
    sample = np.random.Generator(np.random.PCG64(sample_seed))
    labels = [f'S{i+1:02d}' for i in range(P)]
    groups = np.repeat(np.arange(4), [6,6,5,5])
    modules = [f'Constructed module {"ABCD"[i]}' for i in groups]
    # Positive directed ring ensures every node participates; all further edges
    # are sampled without replacement. Orientation is B[outcome, source].
    positive = [((i+1)%P, i) for i in range(P)]
    remaining = [(j,i) for j in range(P) for i in range(P)
                 if j != i and (j,i) not in positive]
    extra = graph.choice(len(remaining), size=31, replace=False)
    positive += [remaining[i] for i in extra[:25]]
    negative = [remaining[i] for i in extra[25:]]
    b = np.diag(np.linspace(.43,.71,P))
    for j,i in positive:
        b[j,i] = graph.uniform(.055,.115)
    for j,i in negative:
        b[j,i] = -graph.uniform(.02,.05)
    # Weak general factor (.15) plus one stronger module factor (.45).
    factors = np.zeros((P,5)); factors[:,0] = .15
    factors[np.arange(P),groups+1] = .45
    correlation = factors@factors.T + np.diag(1-np.sum(factors**2,axis=1))
    sd = np.linspace(.63,.83,P)
    sigma = sd[:,None]*correlation*sd[None,:]
    mu = np.linspace(1.44,2.64,P)
    residual_sd = np.linspace(.38,.52,P)
    psi = np.diag(residual_sd**2)
    intercept = (np.eye(P)-b)@mu
    rho = float(max(abs(np.linalg.eigvals(b))))
    assert rho < 1, 'Constructed system must be stable; do not silently rescale it'
    model = dict(source_labels=labels,outcome_labels=labels,B=b.tolist(),
                 mu1=mu.tolist(),sigma1=sigma.tolist(),psi=psi.tolist(),
                 intercept=intercept.tolist(),anchors=[0]*P,modules=modules,
                 orientation='outcome_by_source',higher_is_worse=True,
                 scale='SYNTHETIC ILLUSTRATIVE raw units; constructed 22-symptom system, NOT original manuscript parameters or clinical evidence')
    scoring = dict(candidates=labels,partners={s:[t for t in labels if t!=s] for s in labels},
                   outcome_weights=[1]*P,utility_weights=[1]*7,horizon=3,gamma=.7,
                   q=.8,breadth_threshold=.05,module_threshold=.03)
    metadata = dict(kind='synthetic_illustrative', original_replication_status='NOT_RUN',
                    seed=SEED,bit_generator='PCG64',child_streams=['graph','paired_sample'],
                    n=N,p=P,module_sizes=[6,6,5,5],
                    ordered_labels=labels,module_map=dict(zip(labels,modules)),
                    symptom_names={s:f'Synthetic symptom {i+1:02d}' for i,s in enumerate(labels)},
                    constructed_choices='All labels, memberships, coefficients, factors, samples, partner sets, H=3, gamma=.7, six doses, and strategy choices are invented for implementation illustration; not recovered original inputs.',
                    spectral_radius=rho,nonzero_offdiagonal=53,positive_offdiagonal=47,negative_offdiagonal=6,
                    general_factor_loading=.15,module_factor_loading=.45,
                    dose_grid=[0,.2,.4,.6,.8,1],
                    interpretation='One n=250 illustrative paired sample; not the manuscript 600-dataset finite-sample experiment, not clinical truth. Spectral radius and ranks were NOT tuned to manuscript results.')
    x = sample.multivariate_normal(mu,sigma,N)
    y = intercept+x@b.T+sample.normal(size=(N,P))*residual_sd
    # Keep paired JSON restricted to the adapter API keys. The adjacent metadata
    # and model scale carry explicit provenance without changing the public API.
    paired = dict(t1=x.tolist(),t2=y.tolist(),
                  fitting=dict(labels=labels,modules=modules,raw_anchors=[0]*P,higher_is_worse=True),scoring=scoring)
    return dict(model=model,scoring=scoring,provenance=metadata,
                dose_requests=[dict(doses={labels[0]:d}) for d in metadata['dose_grid']]),paired,metadata


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent)
    args=parser.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    for name,value in zip(['synthetic_model.json','synthetic_paired.json','synthetic22_metadata.json'],generate()):
        (args.output_dir/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    print('Generated SYNTHETIC ILLUSTRATIVE 22-symptom / 4-module model and 250 paired participants; original replication NOT_RUN')

if __name__=='__main__':main()
