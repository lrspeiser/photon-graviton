"""CC-1 local common-cone Hamiltonian and spherical 3D regulator."""
import numpy as np


def geometry(f, g, eta, kappa):
    alpha = np.exp(g*f[0])
    a = f[1:]
    scale = np.sqrt(1+eta**2*np.dot(a, a))
    b = kappa*eta*a/scale
    db = kappa*eta*(np.eye(3)/scale-eta**2*np.outer(a, a)/scale**3)
    return alpha, alpha*b, db


def particle(f, p, mass, g, eta, kappa):
    alpha, beta, db = geometry(f, g, eta, kappa)
    e = np.sqrt(mass**2+np.dot(p, p))
    h = alpha*e+np.dot(beta, p)
    v = alpha*p/e+beta
    source = np.concatenate(([g*h], alpha*db@p))
    hessian = alpha*(np.eye(3)/e-np.outer(p, p)/e**3)
    return h, v, source, hessian


def spherical_weights(x, q, radius):
    d = x-q
    z = np.maximum(1-np.sum(d*d, axis=-1)/radius**2, 0)
    raw = z**3
    draw = 6*d/radius**2*z[..., None]**2
    norm = raw.sum()
    if norm == 0: raise ValueError('Unresolved source')
    w = raw/norm
    dw = (draw-w[..., None]*draw.sum(axis=(0, 1, 2)))/norm
    return w, dw


def rotation(angle):
    axis = np.array([1., 2., 3.])/np.sqrt(14)
    x, y, z = axis
    cross = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])
    return np.eye(3)*np.cos(angle)+(1-np.cos(angle))*np.outer(axis, axis)+np.sin(angle)*cross
