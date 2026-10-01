"""Illustrate explicit state/transition/strategy distinctions on synthetic data."""
import json
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'skills/tsymperturb/scripts'))
from tsymperturb import (Model, acquisition_order, block, bounded_normal_mean,
                        horizon_utility, reference_block_response, serializable,
                        trajectories)


def main():
    data=json.loads((ROOT/'examples/synthetic_model.json').read_text())
    m=Model(**data['model'])
    delta=np.zeros(m.p);delta[0]=m.mu1[0]-m.anchors[0]
    # Explicit assumption: this illustrative simulation reuses one B, not observed waves.
    path=trajectories([m.B]*3,delta)
    u=np.zeros((3,m.p));u[1,1]=.2
    repeated=trajectories([m.B]*3,delta,interventions=u)
    raw=m.state({'S1':.5})
    bounded_baseline=bounded_normal_mean(m.mu2,m.sd2,0,4)
    bounded_perturbed=bounded_normal_mean(raw['outcome_mean'],np.sqrt(np.diag(raw['outcome_covariance'])),0,4)
    result={
        'stationarity_assumption':'Same B reused beyond one observed interval; strategy simulation, not identified treatment order.',
        'one_time_response_path_raw':path,
        'discounted_raw_utility':horizon_utility(path,[1]*m.p,.9),
        'repeated_additional_improvement':u,
        'repeated_response_path_raw':repeated,
        'signed_edge_block_at_explicit_mean_profile':reference_block_response(m,block(m,'S4',.8,outcome='S1'),m.mu1),
        'bounded_half_dose_response':bounded_baseline-bounded_perturbed,
        'bounded_assumption':'Gaussian marginal clipping to synthetic [0,4]; extra nonlinear operation, not linear reference tVPPS.',
        'target_acquisition':acquisition_order({'S1':.4,'S2':.3,'S3':.2},{'S1':.02,'S2':.01,'S3':.01},2,.9),
    }
    print(json.dumps(serializable(result),indent=2,allow_nan=False))

if __name__=='__main__':main()
