"""FM-1 semidiscrete Hamiltonian; all quantities dimensionless."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


class Model:
    def __init__(self, spec, n=48):
        self.spec, self.n = spec, n
        self.dx, self.a, self.omega = 16/n, .75, .2
        self.g, self.eta = spec['g'], spec['eta']
        x = np.arange(n)*self.dx-8
        self.x, self.y = np.meshgrid(x, x, indexing='ij')
        self.m = np.ones(6 if spec['probe'] == 'none' else 7)
        if spec['probe'] != 'none':
            self.m[-1] = 0 if spec['probe'] == 'photon' else 1e-4
        self.k, self.nf = len(self.m), 3*n*n

    def unpack(self, state):
        f = state[:self.nf].reshape(3, self.n, self.n)
        pi = state[self.nf:2*self.nf].reshape(f.shape)
        q = state[2*self.nf:2*self.nf+2*self.k].reshape(self.k, 2)
        p = state[2*self.nf+2*self.k:].reshape(self.k, 2)
        return f, pi, q, p

    def initial(self):
        theta = np.arange(6)*2*np.pi/6
        q = .8*np.column_stack((np.cos(theta), np.sin(theta)))
        direction = {'rest': 0, 'rotating': 1, 'reverse': -1}[self.spec['source']]
        p = .2*direction*np.column_stack((-np.sin(theta), np.cos(theta)))
        if self.spec['probe'] != 'none':
            photon = self.spec['probe'] == 'photon'
            q = np.vstack((q, [-2.5 if photon else -1.5, .7]))
            momentum = self.spec.get('photon_p', 1e-4) if photon else 1e-4*.4/np.sqrt(1-.4**2)
            p = np.vstack((p, [momentum, 0]))
        return np.concatenate((np.zeros(2*self.nf), q.ravel(), p.ravel()))

    def weights(self, q):
        d = (np.stack((self.x, self.y))[None]-q[:, :, None, None]+8)%16-8
        z = np.maximum(1-(d/self.a)**2, 0)
        b, db = z**3, 6*d/self.a**2*z**2
        raw = b[:, 0]*b[:, 1]
        draw = np.stack((db[:, 0]*b[:, 1], b[:, 0]*db[:, 1]), axis=1)
        norm = raw.sum(axis=(1, 2))
        dnorm = draw.sum(axis=(2, 3))
        w = raw/norm[:, None, None]
        dw = (draw-w[:, None]*dnorm[:, :, None, None])/norm[:, None, None, None]
        return w, dw

    def particle(self, f, q, p):
        w, dw = self.weights(q)
        bar = np.einsum('aij,kij->ka', f, w)
        grad = np.einsum('aij,kbij->kab', f, dw)
        e = np.sqrt(self.m**2+np.sum(p*p, axis=1))
        if np.any(e == 0):
            raise ValueError('Undefined zero-energy massless direction')
        u = p/e[:, None]
        au = np.sum(bar[:, 1:]*u, axis=1)
        exp = np.exp(self.g*bar[:, 0]+self.eta*au)
        h = e*exp
        velocity = exp[:, None]*(u+self.eta*(bar[:, 1:]-u*au[:, None]))
        force = -h[:, None]*(self.g*grad[:, 0]+self.eta*np.einsum('kab,ka->kb', grad[:, 1:], u))
        coeff = np.column_stack((np.full(self.k, self.g), self.eta*u))
        return h, velocity, force, w, coeff

    def lap(self, f):
        return sum(np.roll(f, 1, axis=a)+np.roll(f, -1, axis=a)-2*f for a in (-2, -1))/self.dx**2

    def rhs(self, state):
        f, pi, q, p = self.unpack(state)
        h, v, force, w, coeff = self.particle(f, q, p)
        dpi = self.lap(f)-self.omega**2*f-np.einsum('k,ka,kij->aij', h, coeff, w)/self.dx**2
        out = np.concatenate((pi.ravel(), dpi.ravel(), v.ravel(), force.ravel()))
        if not np.all(np.isfinite(out)):
            raise FloatingPointError('Nonfinite Hamiltonian evolution')
        return out

    def energies(self, state):
        f, pi, q, p = self.unpack(state)
        gradient = sum(((np.roll(f, -1, axis=a)-f)/self.dx)**2 for a in (-2, -1))
        field = .5*self.dx**2*np.sum(pi*pi+gradient+self.omega**2*f*f)
        matter = self.particle(f, q, p)[0].sum()
        return np.array([field+matter, field, matter])

    def diagnostics(self, state):
        f, pi, q, p = self.unpack(state)
        gx = (np.roll(f, -1, axis=1)-np.roll(f, 1, axis=1))/(2*self.dx)
        gy = (np.roll(f, -1, axis=2)-np.roll(f, 1, axis=2))/(2*self.dx)
        momentum = p.sum(axis=0)-self.dx**2*np.array([np.sum(pi*gx), np.sum(pi*gy)])
        angular = np.sum(q[:, 0]*p[:, 1]-q[:, 1]*p[:, 0])
        angular -= self.dx**2*np.sum(pi*(self.x*gy-self.y*gx))
        angular += self.dx**2*np.sum(f[1]*pi[2]-f[2]*pi[1])
        curl = gx[2]-gy[1]
        v = self.particle(f, q, p)[1]
        return dict(momentum=momentum, angular=angular, max_speed=float(np.linalg.norm(v, axis=1).max()),
                    amplitude=float(np.max(abs(f))), curl_max=float(np.max(abs(curl))),
                    circulation_inner=float(self.dx**2*np.sum(curl[self.x**2+self.y**2<2**2])),
                    probe_angle=float(np.arctan2(v[-1, 1], v[-1, 0])) if self.k == 7 else None)


def rk4(rhs, y, dt):
    k1 = rhs(y)
    k2 = rhs(y+.5*dt*k1)
    k3 = rhs(y+.5*dt*k2)
    k4 = rhs(y+dt*k3)
    return y+dt/6*(k1+2*k2+2*k3+k4)


def controls():
    spec = dict(g=.08, eta=.08, source='rotating', probe='photon')
    model = Model(spec, 24)
    y = model.initial()
    rng = np.random.default_rng(20260920)
    y[:2*model.nf] = rng.normal(0, .003, 2*model.nf)
    rhs = model.rhs(y)
    expected = np.concatenate((-rhs[model.nf:2*model.nf]*model.dx**2,
                               rhs[:model.nf]*model.dx**2,
                               -rhs[2*model.nf+2*model.k:], rhs[2*model.nf:2*model.nf+2*model.k]))
    rows = []
    def add(name, error, tolerance):
        rows.append(dict(name=name, error=float(error), tolerance=tolerance, passed=bool(error<=tolerance)))
    # Include every particle coordinate, plus nonzero random field/Pi coordinates.
    indices = list(rng.choice(2*model.nf, 32, replace=False))+list(range(2*model.nf, len(y)))
    for j in indices:
        step = 1e-8 if j >= 2*model.nf+2*model.k and j >= len(y)-2 else 1e-6
        plus, minus = y.copy(), y.copy()
        plus[j] += step
        minus[j] -= step
        numerical = (model.energies(plus)[0]-model.energies(minus)[0])/(2*step)
        add(f'H_gradient_{j}', abs(numerical-expected[j])/max(1, abs(expected[j])), 1e-5)
    f, pi, q, p = model.unpack(y)
    h, v, force, w, coeff = model.particle(f, q, p)
    add('weight_normalization', np.max(abs(w.sum(axis=(1, 2))-1)), 1e-12)
    add('weight_gradient_normalization', np.max(abs(model.weights(q)[1].sum(axis=(2, 3)))), 1e-12)
    dep = np.einsum('k,ka,kij->aij', h, coeff, w)/model.dx**2
    add('source_deposition_integral', np.max(abs(dep.sum(axis=(1, 2))*model.dx**2-(h[:, None]*coeff).sum(axis=0))), 1e-12)
    grad_energy = sum(np.sum(((np.roll(f, -1, axis=a)-f)/model.dx)**2) for a in (-2, -1))
    add('integration_by_parts', abs(np.sum(f*model.lap(f))+grad_energy), 1e-12)
    free = Model(dict(spec, g=0, eta=0), 24)
    r = free.rhs(free.initial())
    add('zero_coupling_field_and_force', max(np.max(abs(r[:2*free.nf])), np.max(abs(r[2*free.nf+2*free.k:]))), 1e-14)
    add('zero_source_zero_field', np.max(abs(model.lap(np.zeros_like(f)))), 1e-14)
    add('positive_field_and_particle_energy', 0 if model.energies(y)[1]>0 and np.all(h>0) else 1, 0)
    return rows


def wave_controls():
    rows, paths = [], {}
    for n in (128, 256, 512):
        dx = 16/n
        x = np.arange(n)*dx-8
        pulse = lambda z: np.maximum(1-(z/.75)**2, 0)**4
        y = np.concatenate((pulse(x), np.zeros(n)))
        dt = 2/np.ceil(2/(.2*dx))
        def rhs(z):
            f, pi = z[:n], z[n:]
            return np.concatenate((pi, (np.roll(f, 1)+np.roll(f, -1)-2*f)/dx**2))
        for _ in range(round(2/dt)):
            y = rk4(rhs, y, dt)
        exact = .5*(pulse(x-2)+pulse(x+2))
        error = np.max(abs(y[:n]-exact))/np.max(exact)
        ahead = np.max(abs(y[:n][abs(x)>.75+2+2*dx]))
        rows.append(dict(n=n, dt=dt, relative_max_error=float(error), ahead_amplitude=float(ahead),
                         error_passed=bool(error<=.03), front_passed=bool(ahead<=.001)))
        paths[str(n)] = np.column_stack((x, y[:n], exact))
    return rows, paths
