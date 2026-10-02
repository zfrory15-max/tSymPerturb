"""Illustrate state/transition/strategy distinctions on a constructed 22-symptom system."""
import json
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'skills/tsymperturb/scripts'))
from tsymperturb import (Model, acquisition_order, block, bounded_normal_mean,
                        horizon_utility, reference_block_response, serializable,
                        trajectories, score)


def main():
    data=json.loads((ROOT/'examples/synthetic_model.json').read_text())
    m=Model(**data['model'])
    delta=np.zeros(m.p);delta[0]=m.mu1[0]-m.anchors[0]
    # Explicit assumption: this illustrative simulation reuses one B, not observed waves.
    path=trajectories([m.B]*4,delta)
    u=np.zeros((4,m.p));u[1,1]=.2
    repeated=trajectories([m.B]*4,delta,interventions=u)
    raw=m.state({'S01':.5})
    bounded_baseline=bounded_normal_mean(m.mu2,m.sd2,0,4)
    bounded_perturbed=bounded_normal_mean(raw['outcome_mean'],np.sqrt(np.diag(raw['outcome_covariance'])),0,4)
    scores=score(m,**data['scoring'])
    top8=np.argsort(-scores['tvpps'],kind='stable')[:8]
    selected=[m.labels[i] for i in top8]
    utilities={m.labels[i]:float(scores['raw'][i,0]) for i in top8}
    costs={name:.01 for name in selected}
    result={
        'provenance':data['provenance'],
        'dimensions':{'symptoms':m.p,'modules':len(set(m.modules)),'transitions':4},
        'acquisition_candidate_subset':selected,
        'acquisition_subset_rule':'Top eight synthetic population tVPPS among all 22, stable label-order tie break; exact two-target search (56 orders), equal constructed cost .01 and discount .9; no clinical interpretation.',
        'stationarity_assumption':'Same B reused beyond one observed interval; strategy simulation, not identified treatment order.',
        'one_time_response_path_raw':path,
        'discounted_raw_utility':horizon_utility(path,[1]*m.p,.9),
        'repeated_additional_improvement':u,
        'repeated_response_path_raw':repeated,
        'signed_edge_block_at_explicit_mean_profile':reference_block_response(m,block(m,'S01',.8,outcome='S02'),m.mu1),
        'bounded_half_dose_response':bounded_baseline-bounded_perturbed,
        'bounded_assumption':'Gaussian marginal clipping to synthetic [0,4]; extra nonlinear operation, not linear reference tVPPS.',
        'target_acquisition':acquisition_order(utilities,costs,2,.9),
    }
    print(json.dumps(serializable(result),indent=2,allow_nan=False))

if __name__=='__main__':main()
