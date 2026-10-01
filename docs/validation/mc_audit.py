import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'skills'/'tsymperturb'/'scripts'))
import tsymperturb as t
from oracle import *
rng=np.random.default_rng(20260930);n=120000
x=rng.multivariate_normal(mu,Sigma,n);noise=rng.multivariate_normal(np.zeros(4),Psi,n)
y=x@B.T+noise
maxz=0
for i in range(4):
 xx=x.copy();xx[:,i]=anchor[i]
 yy=xx@B.T+noise
 diff=y-yy
 observed=diff.mean(axis=0)
 sem=diff.std(axis=0,ddof=1)/np.sqrt(n)
 z=np.max(np.abs((observed-R[:,i])/sem));maxz=max(maxz,z)
 assert z<5,(i,z)
 f=np.ones(4);f[i]=0
 expected=B@np.diag(f)@Sigma@np.diag(f)@B.T+Psi
 assert np.max(np.abs(np.cov(yy,rowvar=False)-expected))<.02
print(f'PASS Monte Carlo: n={n}, seed=20260930, largest mean-error/MCSE={maxz:.6f}; covariance max tolerance .02 raw units squared')
