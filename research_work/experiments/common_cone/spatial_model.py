"""SR-1 spatial response, all sectors derived from one discrete Hamiltonian."""
import numpy as np
from evolution import Evolution
from model import spherical_weights


def coefficients(f, s, lam, g, eta, kappa):
    u = np.concatenate((np.full_like(f[:1], g), -lam*f[1:]), axis=0)
    potential = g*f[0]-lam*np.sum(f[1:]**2, axis=0)/2
    alpha = np.exp(potential); z = np.exp(2*s*potential)
    a = np.exp((1+3*s)*potential); d = np.exp((1-s)*potential)
    c = np.exp((1+s)*potential)
    scale = np.sqrt(1+eta**2*np.sum(f[1:]**2, axis=0))
    b = kappa*eta*f[1:]/scale; beta = c*b
    shape = (3, 4)+f.shape[1:]
    db = np.zeros(shape)
    for j in range(3):
        for k in range(3):
            db[j,k+1] = kappa*eta*((j == k)/scale-eta**2*f[j+1]*f[k+1]/scale**3)
    B = (1+s)*beta[:,None]*u[None]+c*db
    return u, alpha, z, a, d, c, beta, B


class SpatialEvolution(Evolution):
    def __init__(self, s=1., lam=0., source_momentum=.2, **kwargs):
        super().__init__(**kwargs)
        self.s, self.lam, self.source_momentum = s, lam, source_momentum

    def initial(self):
        y = super().initial()
        _,_,_,p = self.unpack(y)
        p[:6] *= self.source_momentum/.2
        return y

    def pieces(self,y):
        f,pi,q,p = self.unpack(y)
        u,alpha,z,a,d,c,beta,B = coefficients(f,self.s,self.lam,self.g,self.eta,self.kappa)
        gp = np.stack([self.plus(f,j+1) for j in range(3)])
        gm = np.stack([self.minus(f,j+1) for j in range(3)])
        gc = (gp+gm)/2
        kinetic = np.sum(pi*pi,axis=0)/2
        gradient = np.sum(gp*gp+gm*gm,axis=(0,1))/4
        current = np.sum(pi[None]*gc,axis=1)
        field_density = a*kinetic+d*gradient-np.sum(beta*current,axis=0)+self.omega**2*np.sum(f*f,axis=0)/2
        w,dw = zip(*(spherical_weights(self.x,pos,self.radius) for pos in q))
        w,dw = np.array(w),np.array(dw)
        p2 = np.sum(p*p,axis=1)[:,None,None,None]
        e = np.sqrt(self.m[:,None,None,None]**2+z*p2)
        local_h = alpha*e+np.einsum('jxyt,kj->kxyt',beta,p)
        velocity = p*np.sum(w*alpha*z/e,axis=(1,2,3))[:,None]+np.einsum('kxyz,jxyz->kj',w,beta)
        force = -np.einsum('kxyzj,kxyz->kj',dw,local_h)
        return locals()

    def rhs(self,y,omit_self=False):
        t = self.pieces(y)
        f,pi,u,a,d,beta,B,gp,gm,gc = (t[k] for k in ('f','pi','u','a','d','beta','B','gp','gm','gc'))
        fd = a*pi-np.sum(beta[:,None]*gc,axis=0)
        pd = -self.omega**2*f
        for axis in range(3):
            pd += .5*(self.minus(d*gp[axis],axis+1)+self.plus(d*gm[axis],axis+1))
            pd -= self.central(beta[axis]*pi,axis+1)
        if not omit_self:
            pd -= u*((1+3*self.s)*a*t['kinetic']+(1-self.s)*d*t['gradient'])
            pd += np.einsum('jaxyz,jxyz->axyz',B,t['current'])
        scalar = np.sum(t['w']*t['alpha']*(t['e']+self.s*t['z']*t['p2']/t['e']),axis=0)/self.dv
        momentum = np.einsum('kxyz,kj->jxyz',t['w'],t['p'])/self.dv
        pd -= u*scalar+np.einsum('jaxyz,jxyz->axyz',B,momentum)
        pd -= self.gamma*fd
        loss = np.sum(self.gamma*fd*fd)*self.dv
        out = np.concatenate((fd.ravel(),pd.ravel(),t['velocity'].ravel(),t['force'].ravel(),[loss]))
        if not np.all(np.isfinite(out)): raise FloatingPointError('SR-1 derivative nonfinite')
        return out

    def metrics(self,y):
        t = self.pieces(y)
        field = np.sum(t['field_density'])*self.dv
        matter = np.sum(t['w']*t['local_h'])
        bbar = np.einsum('kxyz,jxyz->kj',t['w'],t['beta'])
        cbar = np.sum(t['w']*t['c'],axis=(1,2,3))
        relative = np.linalg.norm(t['velocity']-bbar,axis=1)/cbar
        cone = max(0.,float(relative.max()-1))
        if self.probe == 'photon': cone = max(cone,abs(float(relative[-1]-1)))
        j = t['current']; f = t['f']; pi = t['pi']; p = t['p']; q = t['q']
        momentum = p.sum(axis=0)-np.sum(j,axis=(1,2,3))*self.dv
        angular = np.cross(q,p).sum(axis=0)-np.cross(self.x,j.transpose(1,2,3,0)).sum(axis=(0,1,2))*self.dv
        angular += np.cross(f[1:].transpose(1,2,3,0),pi[1:].transpose(1,2,3,0)).sum(axis=(0,1,2))*self.dv
        curlz = t['gc'][0,2]-t['gc'][1,1]
        plane = self.n//2; mask = np.sum(self.x[:,:,plane,:2]**2,axis=-1)<4
        edge = max(np.max(abs(np.take(f,i,axis=axis))) for axis in (1,2,3) for i in (0,-1))
        return dict(total=field+matter,field=field,matter=matter,absorbed=float(y[-1]),ledger=field+matter+y[-1],
                    momentum=momentum,angular=angular,cone_residual=cone,
                    max_characteristic=float(np.max(t['c']+np.linalg.norm(t['beta'],axis=0))),
                    edge_amplitude=float(edge),source_extent=float(np.max(abs(q))+self.radius),
                    circulation_z=float(np.sum(curlz[:,:,plane][mask])*self.dx**2),
                    probe_angle=float(np.arctan2(t['velocity'][-1,1],t['velocity'][-1,0])))


def local(f,pi,grad,p,m,s,lam,g=.08,eta=.08,kappa=.5,omega=.2):
    u,alpha,z,a,d,c,beta,B = coefficients(f,s,lam,g,eta,kappa)
    e = np.sqrt(m*m+z*np.dot(p,p)); j = grad@pi
    hf = a*np.dot(pi,pi)/2+d*np.sum(grad*grad)/2-np.dot(beta,j)+omega**2*np.dot(f,f)/2
    hp = alpha*e+np.dot(beta,p)
    df = u*((1+3*s)*a*np.dot(pi,pi)/2+(1-s)*d*np.sum(grad*grad)/2)+omega**2*f-B.T@j
    df += u*alpha*(e+s*z*np.dot(p,p)/e)+B.T@p
    dpi = a*pi-beta@grad
    dg = d*grad-beta[:,None]*pi
    dp = alpha*z*p/e+beta
    hpp = alpha*z*(np.eye(3)/e-z*np.outer(p,p)/e**3)
    return hf+hp,df,dpi,dg,dp,hpp,(a,d,c,beta)
