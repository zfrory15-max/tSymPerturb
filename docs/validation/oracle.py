"""Independent manuscript-equation oracle; does not call implementation."""
import numpy as np
B=np.array([[.5,.2,-.1,.04],[.3,.4,.05,-.02],[-.1,.1,.6,.2],[.08,-.03,.25,.45]])
mu=np.array([2.,1.,3.,2.5]); anchor=np.array([0.,-.5,1.,0.])
Sigma=np.array([[1.,.2,.1,0.],[.2,2.,.3,.1],[.1,.3,1.5,.2],[0.,.1,.2,.8]])
Psi=np.diag([.4,.3,.2,.5]); w=np.array([1.,2.,3.,4.]); modules=np.array(['A','A','B','B'])
sd=np.sqrt(np.diag(B@Sigma@B.T+Psi)); R=B*(mu-anchor)[None,:]; Z=R/sd[:,None]
def q(mat,H=3,gamma=.7):
 return sum(np.linalg.matrix_power(gamma*np.abs(mat),h).sum() for h in range(1,H+1))
raw=[]
for i in range(4):
 other=np.arange(4)!=i
 cb=B.copy(); cb[:,i]*=.2
 pair=[]
 for k in range(4):
  if k==i: continue
  rest=(np.arange(4)!=i)&(np.arange(4)!=k)
  gi=np.average(Z[rest,i],weights=w[rest]); gk=np.average(Z[rest,k],weights=w[rest]); gij=np.average(Z[rest,i]+Z[rest,k],weights=w[rest]); pair.append(gij-max(gi,gk))
 raw.append([np.average(Z[:,i],weights=w),np.average(Z[other,i],weights=w[other]),np.mean(Z[other,i]>=.05),np.mean([Z[modules==m,i].mean()>=.03 for m in set(modules) if m!=modules[i]]),(q(B)-q(cb))/q(B),np.mean(np.maximum(pair,0)),sum(w[other]*np.maximum(Z[other,i],0))/(sum(w*np.maximum(Z[:,i],0))+1e-12)])
raw=np.array(raw); spans=np.ptp(raw,axis=0); scores=np.empty_like(raw)
for j in range(7): scores[:,j]=50 if spans[j]==0 else 100*(raw[:,j]-raw[:,j].min())/spans[j]
if __name__=='__main__':
 print('sd',sd);print('R',R);print('raw seven',raw);print('score',scores.mean(axis=1))

