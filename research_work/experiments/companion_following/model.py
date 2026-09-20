"""GF-1 finite packet mechanics; see protocol.md for scope and attribution."""
from __future__ import annotations
import hashlib
import itertools
import json
import numpy as np

TAU = 0.5
KERNELS = ["exp(-s^2)", "exp(-s)", "(1+s^2)^(-1)", "(1+s^2)^(-2)",
           "max(1-s,0)^4*(1+4*s)", "s^2*exp(1-s^2)"]
GATES = ["1", "(1+c)/2", "max(c,0)", "max(c,0)^2", "max(c,0)*(1+h)/2"]
MEMORIES = ["u=v", "tau*dz1/dt=v-z1; u=z1",
            "tau*dz1/dt=v-z1; tau*dz2/dt=z1-z2; u=z2",
            "tau^2*d2z2/dt2+0.7*tau*dz2/dt+z2=v; u=z2"]
FORCES = [
    "sum_j W_ij*(u_j-v_i)",
    "P_i*sum_j W_ij*(u_j-v_i)",
    "P_i*sum_j W_ij*(u_j-v_i+tau*omega_j*J*u_j)",
    "sum_j S_ij*cross(u_i-u_j,e_ij)*J*(v_i-v_j)",
    "sum_j W_ij*(d_ij-tau*u_j)/(1+|d_ij-tau*u_j|^2)^(3/2)",
]

def catalogue():
    rows = []
    for f, k, g, m in itertools.product(range(5), range(6), range(5), range(4)):
        equation = (f"a_i=-grad_i U+0.2*({FORCES[f]}); "
                    f"W_ij=({KERNELS[k]})*({GATES[g]}); {MEMORIES[m]}")
        if f == 2:
            equation += "; tau*domega_i/dt=cross(v_i,a_i)/|v_i|^2-omega_i"
        rows.append(dict(id=f"F{f}-K{k}-G{g}-M{m}", family=f, kernel=k,
                         gate=g, memory=m, eta=0., equation=equation,
                         status="phenomenological; instantaneous spatial interaction"))
    for k, eta in itertools.product(range(6), (.5, 2., 8.)):
        equation = (f"L=0.5*sum_i|v_i|^2+{eta}/32*sum_(i<j)"
                    f"({KERNELS[k]})*|v_i-v_j|^2-U; M*a=b (protocol Euler-Lagrange)")
        rows.append(dict(id=f"C-K{k}-E{eta:g}", family=5, kernel=k,
                         gate=0, memory=0, eta=eta, equation=equation,
                         status="conservative finite packet Lagrangian; instantaneous interaction"))
    for r in rows:
        r["equation_sha256"] = hashlib.sha256(r["equation"].encode()).hexdigest()
    return rows

def cross(a, b):
    return a[..., 0]*b[..., 1]-a[..., 1]*b[..., 0]

def J(a):
    return np.stack((-a[..., 1], a[..., 0]), axis=-1)

def unit(a):
    norm = np.linalg.norm(a, axis=-1, keepdims=True)
    return np.divide(a, norm, out=np.zeros_like(a), where=norm != 0)

def geometry(x):
    d = x[:, None, :, :]-x[:, :, None, :]
    r = np.linalg.norm(d, axis=-1)
    e = np.divide(d, r[..., None], out=np.zeros_like(d), where=r[..., None] != 0)
    mask = 1-np.eye(x.shape[1])[None, :, :]
    return d, r, e, mask

def kernel(r, ids):
    s = r/2.5
    out = np.empty_like(s)
    deriv = np.empty_like(s)
    for i in range(6):
        ix = ids == i
        z = s[ix]
        if not np.any(ix):
            continue
        if i == 0:
            a = np.exp(-z*z); b = -2*z*a
        elif i == 1:
            a = np.exp(-z); b = -a
        elif i == 2:
            a = 1/(1+z*z); b = -2*z*a*a
        elif i == 3:
            a = (1+z*z)**-2; b = -4*z*(1+z*z)**-3
        elif i == 4:
            q = np.maximum(1-z, 0)
            a = q**4*(1+4*z); b = -20*z*q**3
        else:
            a = z*z*np.exp(1-z*z)
            b = 2*z*(1-z*z)*np.exp(1-z*z)
        out[ix] = a
        deriv[ix] = b/2.5
    return out, deriv

def potential_force(r, e, mask):
    short = np.exp(-r/.4)
    long = np.exp(-r/2)
    U = np.sum((2*short-long)*mask, axis=(1, 2))/32
    force = np.sum(((-5*short+.5*long)*mask)[..., None]*e, axis=2)/16
    return U, force

class Mechanics:
    def __init__(self, rows, coupling=.2):
        self.rows = rows
        self.f = np.array([r["family"] for r in rows])
        self.k = np.array([r["kernel"] for r in rows])
        self.g = np.array([r["gate"] for r in rows])
        self.m = np.array([r["memory"] for r in rows])
        self.eta = np.array([r["eta"] for r in rows])
        self.conservative = bool(np.all(self.f == 5))
        if not (self.conservative or np.all(self.f < 5)):
            raise ValueError("Use separate conservative/phenomenological batches")
        self.coupling = coupling

    def terms(self, state):
        x, v = state[..., :2], state[..., 2:4]
        d, r, e, mask = geometry(x)
        U, fc = potential_force(r, e, mask)
        K, kr = kernel(r, self.k)
        K *= mask
        kr *= mask
        if self.conservative:
            n = x.shape[1]
            lap = -K.copy()
            ii = np.arange(n)
            lap[:, ii, ii] += np.sum(K, axis=2)
            kk = self.eta[:, None, None]/16
            M = np.eye(n)[None] + kk*lap
            dv = v[:, :, None, :]-v[:, None, :, :]
            er = -e
            bkin = (self.eta[:, None, None]/16)*np.sum(
                kr[..., None]*(.5*er*np.sum(dv*dv, axis=-1)[..., None]
                  -np.sum(er*dv, axis=-1)[..., None]*dv), axis=2)
            p = np.einsum("bij,bjk->bik", M, v)
            energy = .5*np.sum(v*p, axis=(1, 2))+U
            return dict(U=U, fc=fc, M=M, b=fc+bkin, p=p, energy=energy,
                        r=r, guide=bkin, guide_work=np.sum(v*bkin, axis=(1, 2)))
        u = v.copy()
        u[self.m == 1] = state[self.m == 1, :, 4:6]
        u[self.m >= 2] = state[self.m >= 2, :, 6:8]
        n = unit(v)
        c = np.sum(n[:, :, None, :]*e, axis=-1)
        h = np.sum(n[:, :, None, :]*n[:, None, :, :], axis=-1)
        G = np.ones_like(r)
        for g in range(1, 5):
            ix = self.g == g
            if g == 1: G[ix] = (1+c[ix])/2
            if g == 2: G[ix] = np.maximum(c[ix], 0)
            if g == 3: G[ix] = np.maximum(c[ix], 0)**2
            if g == 4: G[ix] = np.maximum(c[ix], 0)*(1+h[ix])/2
        W = K*G
        align = np.einsum("bij,bjk->bik", W, u)-np.sum(W, axis=2)[..., None]*v
        phi = align.copy()
        ix = self.f == 2
        if np.any(ix):
            source = state[ix, :, 8, None]*J(u[ix])
            phi[ix] += TAU*np.einsum("bij,bjk->bik", W[ix], source)
        ix = (self.f == 1) | (self.f == 2)
        phi[ix] -= n[ix]*np.sum(n[ix]*phi[ix], axis=-1)[..., None]
        phi[ix] *= (np.linalg.norm(v[ix], axis=-1) != 0)[..., None]
        ix = self.f == 3
        if np.any(ix):
            du = u[ix, :, None, :]-u[ix, None, :, :]
            dv = v[ix, :, None, :]-v[ix, None, :, :]
            S = (W[ix]+W[ix].transpose(0, 2, 1))/2
            phi[ix] = np.sum((S*cross(du, e[ix]))[..., None]*J(dv), axis=2)
        ix = self.f == 4
        if np.any(ix):
            q = d[ix]-TAU*u[ix, None, :, :]
            phi[ix] = np.sum(W[ix, ..., None]*q/
                (1+np.sum(q*q, axis=-1))[..., None]**1.5, axis=2)
        guide = self.coupling*phi
        return dict(U=U, fc=fc, b=fc+guide, p=v,
                    energy=.5*np.sum(v*v, axis=(1, 2))+U, r=r, guide=guide,
                    guide_work=np.sum(v*guide, axis=(1, 2)))

    def rhs(self, t, state, chain=False, amplitude=.4):
        v = state[..., 2:4]
        z = self.terms(state)
        leader_a = leader_acceleration(t, amplitude)
        if self.conservative:
            if chain:
                a = np.zeros_like(v)
                a[:, 0] = leader_a
                a[:, 1:] = np.linalg.solve(z["M"][:, 1:, 1:],
                    z["b"][:, 1:]-z["M"][:, 1:, 0, None]*leader_a)
                reaction = np.einsum("bi,bij->bj", z["M"][:, 0, :], a)-z["b"][:, 0]
            else:
                a = np.linalg.solve(z["M"], z["b"])
                reaction = np.zeros((len(state), 2))
        else:
            a = z["b"].copy()
            reaction = np.zeros((len(state), 2))
            if chain:
                reaction = leader_a-a[:, 0]
                a[:, 0] = leader_a
        out = np.zeros_like(state)
        out[..., :2] = v
        out[..., 2:4] = a
        ix = (self.m == 1) | (self.m == 2)
        out[ix, :, 4:6] = (v[ix]-state[ix, :, 4:6])/TAU
        ix = self.m == 2
        out[ix, :, 6:8] = (state[ix, :, 4:6]-state[ix, :, 6:8])/TAU
        ix = self.m == 3
        # z1 stores dz2/dt only for M3, initialized at zero.
        out[ix, :, 6:8] = state[ix, :, 4:6]
        out[ix, :, 4:6] = (v[ix]-state[ix, :, 6:8]
                            -2*.35*TAU*state[ix, :, 4:6])/TAU**2
        ix = self.f == 2
        speed2 = np.sum(v[ix]*v[ix], axis=-1)
        rate = np.divide(cross(v[ix], a[ix]), speed2,
                         out=np.zeros_like(speed2), where=speed2 != 0)
        out[ix, :, 8] = (rate-state[ix, :, 8])/TAU
        return out, reaction

def leader_angle(t, amplitude=.4):
    return amplitude*np.sin(np.pi*(t-3)/8)**2 if 3 < t < 11 else 0.

def leader_acceleration(t, amplitude=.4):
    angle = leader_angle(t, amplitude)
    rate = amplitude*np.pi/8*np.sin(2*np.pi*(t-3)/8) if 3 < t < 11 else 0.
    return rate*np.array([-np.sin(angle), np.cos(angle)])

def rk4(mech, t, y, dt, chain=False, amplitude=.4):
    k1 = mech.rhs(t, y, chain, amplitude)[0]
    k2 = mech.rhs(t+dt/2, y+dt*k1/2, chain, amplitude)[0]
    k3 = mech.rhs(t+dt/2, y+dt*k2/2, chain, amplitude)[0]
    k4 = mech.rhs(t+dt, y+dt*k3, chain, amplitude)[0]
    return y+dt*(k1+2*k2+2*k3+k4)/6

def initial(mech, fixture, n=None):
    n = n or (12 if fixture == "chain" else 16)
    b = len(mech.rows)
    state = np.zeros((b, n, 9))
    periods = np.full(b, np.nan)
    if fixture == "chain":
        state[:, :, 0] = -np.arange(n)
        state[:, :, 2] = 1.
    else:
        theta = 2*np.pi*np.arange(n)/n
        radial = np.stack((np.cos(theta), np.sin(theta)), axis=-1)
        tangent = J(radial)
        state[:, :, :2] = 3*radial
        zero = mech.terms(state)
        a0 = np.linalg.solve(zero["M"], zero["b"]) if mech.conservative else zero["fc"]
        f0 = np.mean(np.sum(a0*radial, axis=-1), axis=1)
        if mech.conservative:
            state[:, :, 2:4] = tangent
            one = mech.terms(state)
            a1 = np.linalg.solve(one["M"], one["b"])-a0
            f1 = np.mean(np.sum(a1*radial, axis=-1), axis=1)
        else:
            f1 = np.zeros(b)
        speed2 = -f0/(1/3+f1)
        if np.any(speed2 <= 0):
            raise ValueError("No positive circular balance for "+str(
                [mech.rows[i]["id"] for i in np.flatnonzero(speed2 <= 0)]))
        speed = np.sqrt(speed2)
        state[:, :, :2] = (3+.02*np.cos(2*theta))[None, :, None]*radial
        state[:, :, 2:4] = speed[:, None, None]*tangent
        state[:, :, 8] = speed[:, None]/3
        periods = 2*np.pi*3/speed
    state[:, :, 4:6] = state[:, :, 2:4]
    state[:, :, 6:8] = state[:, :, 2:4]
    state[mech.m == 3, :, 4:6] = 0
    return state, periods
