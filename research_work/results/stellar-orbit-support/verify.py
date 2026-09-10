"""Independent checks of the mixture optimizer and Jacobi derivatives."""
import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm
from run import fit_weights,invariants,Field,DiskGrid,save

rng=np.random.default_rng(9517)
centres=np.array([-3000.,0.,3000.]);truth=np.array([.2,.5,.3])
labels=rng.choice(3,10000,p=truth);sigma=np.full(len(labels),1000.)
y=centres[labels]+rng.normal(0,sigma)
fit=fit_weights(y,sigma,centres)
K=norm.pdf(y[:,None],loc=centres,scale=sigma[:,None])
def fun(w):
    mix=K@w
    return -np.mean(np.log(mix)), -np.mean(K/mix[:,None],axis=0)
reference=minimize(fun,np.ones(3)/3,jac=True,method='SLSQP',bounds=[(1e-8,1)]*3,
                   constraints=[dict(type='eq',fun=lambda w:w.sum()-1,jac=lambda w:np.ones(3))],
                   options=dict(ftol=1e-12,maxiter=500))
assert reference.success and fit['converged']
dw=float(np.max(abs(np.array(fit['weights'])-reference.x)))
assert dw<1e-4
assert np.max(abs(np.array(fit['weights'])-truth))<.02
f=Field(DiskGrid(256,129,None))
points=np.array([[1.,.1,.3,30.,120.,20.],[3.,.7,-.6,-40.,170.,-30.],[8.,.2,1.1,10.,240.,5.]])
J,grad=invariants(f,points)
numeric=[]
for i,h in enumerate([1e-5,1e-5,1e-5,.001,.001,.001]):
    a=points.copy();b=points.copy();a[:,i]+=h;b[:,i]-=h
    numeric.append((invariants(f,a)[0]-invariants(f,b)[0])/(2*h))
error=float(np.max(abs(np.array(numeric).T-grad)/np.maximum(abs(grad),1)))
assert error<1e-5
save('verification.json',dict(EM_weights=fit['weights'],independent_SLSQP_weights=reference.x.tolist(),
    maximum_weight_difference=dw,known_generating_weights=truth.tolist(),
    maximum_Jacobi_gradient_scaled_difference=error,
    inference_scope='Optimizer and derivative checks only; no claim of exact measurement likelihood or adequate orbit coverage'))
print('Mixture and Jacobi derivative checks passed:',dw,error)
