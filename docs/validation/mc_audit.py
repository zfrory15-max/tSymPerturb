"""22-node Monte Carlo implementation check, not original 22x250000 benchmark.

Uses 120000 common-noise paired baseline/post-KO draws per target. That variance
reduction design is deliberate and differs from the unavailable original design.
"""
import json
import sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'skills/tsymperturb/scripts'))
from tsymperturb import Model
spec=json.loads((ROOT/'examples/synthetic_model.json').read_text())['model']
m=Model(**spec); p=m.p
rng=np.random.default_rng(20261002);n=120000
x=rng.multivariate_normal(m.mu1,m.sigma1,n)
noise=rng.multivariate_normal(np.zeros(p),m.psi,n)
y=m.intercept+x@m.B.T+noise
maxz=0.;maxcov=0.
for i in range(p):
 xx=x.copy();xx[:,i]=m.anchors[i]
 yy=m.intercept+xx@m.B.T+noise
 diff=y-yy;observed=diff.mean(axis=0)
 expected=m.B[:,i]*(m.mu1[i]-m.anchors[i])
 sem=diff.std(axis=0,ddof=1)/np.sqrt(n)
 active=np.abs(m.B[:,i])>0
 z=np.max(np.abs((observed[active]-expected[active])/sem[active]));maxz=max(maxz,z)
 assert z<5,(i,z)
 assert np.max(abs(observed[~active]-expected[~active]),initial=0)<1e-12
 f=np.ones(p);f[i]=0
 cov=m.B@(f[:,None]*m.sigma1*f[None,:])@m.B.T+m.psi
 error=np.max(np.abs(np.cov(yy,rowvar=False)-cov));maxcov=max(maxcov,error)
 assert error<.02,(i,error)
print(f'PASS SYNTHETIC22 Monte Carlo: targets={p}, n={n}, seed=20261002, largest mean-error/MCSE={maxz:.6f}; covariance max error={maxcov:.6f}; original replication NOT_RUN')
