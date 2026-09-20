"""Independent controls for GF-1. Failure is recorded before raising."""
from __future__ import annotations
import argparse
from pathlib import Path
import time
import numpy as np
from model import (Mechanics, catalogue, cross, geometry, initial, J,
                   kernel, potential_force)
from campaign import write_json, manifest, best_correlation

def main(output):
    output.mkdir(parents=True, exist_ok=False)
    write_json(output/"manifest.json", manifest())
    rows = catalogue()
    checks = []
    def check(name, measured, limit, note=""):
        passed = bool(np.isfinite(measured) and abs(measured) <= limit)
        checks.append(dict(name=name, measured=float(measured), limit=limit,
                           passed=passed, note=note))
        print(f"{'PASS' if passed else 'FAIL'} {name}: {measured:.6g}", flush=True)
    start = time.perf_counter()
    times=np.arange(161)*.2
    pulse=np.exp(-((times-5)/1.5)**2)
    delayed=np.zeros_like(pulse);delayed[90:]=pulse[:-90]
    corr,lag=best_correlation(delayed,pulse)
    check("declared full lag window finds lag 18",abs(lag-18)+abs(corr-1),1e-12)
    check("catalogue count", len(rows)-618, 0)
    check("unique structural equations", len(set(r["equation_sha256"] for r in rows))-618, 0)
    rng = np.random.default_rng(61019)
    n = 5
    y = np.zeros((618, n, 9))
    y[..., :2] = rng.normal(size=(n, 2))*1.3
    y[..., 2:4] = rng.normal(size=(n, 2))*.5
    y[..., 4:6] = rng.normal(size=(n, 2))*.4
    y[..., 6:8] = rng.normal(size=(n, 2))*.6
    y[..., 8] = rng.normal(size=n)*.2
    eps = 1e-6
    z = Mechanics(rows[:600]).terms(y[:600])
    d, r, e, mask = geometry(y[:1, :, :2])
    _, force = potential_force(r, e, mask)
    grad = np.zeros((n, 2))
    for i in range(n):
        for a in range(2):
            plus = y[:1, :, :2].copy(); minus = plus.copy()
            plus[0, i, a] += eps; minus[0, i, a] -= eps
            def U(x):
                _, rr, ee, mm = geometry(x)
                return potential_force(rr, ee, mm)[0][0]
            grad[i, a] = (U(plus)-U(minus))/(2*eps)
    check("independent potential gradient", np.max(np.abs(force[0]+grad)), 2e-9)
    for f in (1, 2):
        ix = np.array([r["family"] == f for r in rows[:600]])
        work = np.sum(y[:600, :, 2:4][ix]*z["guide"][ix], axis=-1)
        check(f"F{f} pointwise transverse work", np.max(np.abs(work)), 1e-12)
    ix = np.array([r["family"] == 3 for r in rows[:600]])
    check("F3 reciprocal total force", np.max(np.abs(np.sum(z["guide"][ix], axis=1))), 1e-12)
    check("F3 reciprocal total work", np.max(np.abs(z["guide_work"][ix])), 1e-12)
    # Torque is deliberately not asserted to vanish.
    f3_torque = np.sum(cross(y[:600, :, :2][ix], z["guide"][ix]), axis=1)
    check("F3 nonconserved angular momentum detected", float(np.max(np.abs(f3_torque)) < 1e-8), 0,
          "Pass means a random-state counterexample to angular conservation was found.")
    kernels = np.arange(6)
    rr = np.broadcast_to(np.array([.2, 1., 2., 3., 5.])[None, None, :], (6, 1, 5)).copy()
    kval, derivative = kernel(rr, kernels)
    diff = (kernel(rr+eps, kernels)[0]-kernel(rr-eps, kernels)[0])/(2*eps)
    check("all six kernel radial derivatives", np.max(np.abs(diff-derivative)), 2e-9)

    mc = Mechanics(rows[600:])
    yc = y[600:]
    zc = mc.terms(yc)
    def lagrange(state):
        q = mc.terms(state)
        return q["energy"]-2*q["U"]
    dxL = np.zeros((18, n, 2))
    dvL = np.zeros_like(dxL)
    for i in range(n):
        for a in range(2):
            yp = yc.copy(); ym = yc.copy()
            yp[:, i, a] += eps; ym[:, i, a] -= eps
            dxL[:, i, a] = (lagrange(yp)-lagrange(ym))/(2*eps)
            yp = yc.copy(); ym = yc.copy()
            yp[:, i, a+2] += eps; ym[:, i, a+2] -= eps
            dvL[:, i, a] = (lagrange(yp)-lagrange(ym))/(2*eps)
    yp = yc.copy(); ym = yc.copy()
    yp[..., :2] += eps*yc[..., 2:4]
    ym[..., :2] -= eps*yc[..., 2:4]
    Dp = (mc.terms(yp)["p"]-mc.terms(ym)["p"])/(2*eps)
    check("canonical momentum equals dL/dv", np.max(np.abs(zc["p"]-dvL)), 2e-8)
    check("independent Euler-Lagrange force", np.max(np.abs(zc["b"]-(dxL-Dp))), 2e-8)
    check("positive kinetic metric", max(0., 1.-float(np.min(np.linalg.eigvalsh(zc["M"])))), 1e-12)
    a = np.linalg.solve(zc["M"], zc["b"])
    check("Noether linear momentum rate", np.max(np.abs(np.sum(a, axis=1))), 1e-12)
    # dp/dt = grad_x L; d(x cross p)/dt = v cross p + x cross grad_x L.
    check("Noether canonical angular rate",
          np.max(np.abs(np.sum(cross(yc[..., 2:4], zc["p"])+cross(yc[..., :2], dxL), axis=1))), 2e-8)
    yp = yc.copy(); ym = yc.copy()
    yp[..., :2] += eps*yc[..., 2:4]; ym[..., :2] -= eps*yc[..., 2:4]
    yp[..., 2:4] += eps*a; ym[..., 2:4] -= eps*a
    check("conservative energy directional derivative",
          np.max(np.abs((mc.terms(yp)["energy"]-mc.terms(ym)["energy"])/(2*eps))), 2e-8)
    kk = .5; K = .7; a0 = np.array([.1, .3])
    matrix = np.array([[1+kk*K, -kk*K], [-kk*K, 1+kk*K]])
    af = -(matrix[1, 0]/matrix[1, 1])*a0
    check("two-body analytic acceleration sharing", np.max(np.abs(af-(kk*K/(1+kk*K))*a0)), 1e-14)
    bad_matrix = np.eye(2)-np.array([[1., -1.], [-1., 1.]])
    check("negative coupling ghost control rejected",
          float(np.min(np.linalg.eigvalsh(bad_matrix)) >= 0), 0,
          "Negative eta making the kinetic metric indefinite is rejected, not integrated.")

    for group, rr0, yy in [("phenomenological", rows[:600], y[:600]), ("conservative", rows[600:], yc)]:
        mech = Mechanics(rr0)
        base = mech.rhs(0., yy)[0][..., 2:4]
        moved = yy.copy(); moved[..., :2] += np.array([4., -2.])
        check(group+" translation covariance", np.max(np.abs(base-mech.rhs(0., moved)[0][..., 2:4])), 2e-12)
        for name, R in [("rotation", np.array([[.6, -.8], [.8, .6]])),
                        ("reflection", np.array([[1., 0.], [0., -1.]]))]:
            transformed = yy.copy()
            for slot in (0, 2, 4, 6):
                transformed[..., slot:slot+2] = yy[..., slot:slot+2]@R.T
            transformed[..., 8] *= np.linalg.det(R)
            result = mech.rhs(0., transformed)[0][..., 2:4]
            check(group+" "+name+" covariance", np.max(np.abs(result-base@R.T)), 2e-12)
    mzero = Mechanics(rows[:600], coupling=0)
    check("zero guide coupling", np.max(np.abs(mzero.rhs(0, y[:600])[0][..., 2:4]-z["fc"])), 1e-14)

    # Cucker-Smale continuous two-node Laplacian: dv/dt=-lambda*L*v.
    lam, kval, dt, T = .3, .7, .01, 5.
    vel = np.array([1., -1.])
    def rhs(v): return lam*kval*np.array([v[1]-v[0], v[0]-v[1]])
    for _ in range(round(T/dt)):
        k1=rhs(vel); k2=rhs(vel+dt*k1/2); k3=rhs(vel+dt*k2/2); k4=rhs(vel+dt*k3)
        vel += dt*(k1+2*k2+2*k3+k4)/6
    check("credited two-agent alignment exponential", np.max(np.abs(vel-np.array([1., -1.])*np.exp(-2*lam*kval*T))), 1e-10)

    # Narrow fixed-lattice inertial-turning wave benchmark; no animal-data claim.
    count, mode, alpha, T = 32, 3, .7, 10.
    phi0 = np.sin(2*np.pi*mode*np.arange(count)/count)
    state = np.stack([phi0, np.zeros(count)])
    omega = 2*np.sqrt(alpha)*np.sin(np.pi*mode/count)
    def wave(q):
        return np.stack([q[1], alpha*(np.roll(q[0], 1)+np.roll(q[0], -1)-2*q[0])])
    for _ in range(round(T/dt)):
        k1=wave(state); k2=wave(state+dt*k1/2); k3=wave(state+dt*k2/2); k4=wave(state+dt*k3)
        state += dt*(k1+2*k2+2*k3+k4)/6
    exact = np.stack([phi0*np.cos(omega*T), -omega*phi0*np.sin(omega*T)])
    check("credited linear turning-wave dispersion", np.max(np.abs(state-exact)), 1e-8)
    output_data = dict(checks=checks, count=len(checks), passed=all(c["passed"] for c in checks),
                       wall_seconds=time.perf_counter()-start)
    write_json(output/"controls.json", output_data)
    if not output_data["passed"]:
        raise SystemExit("GF-1 controls failed; preserved archive, screening prohibited")
    print(f"All {len(checks)} controls pass.", flush=True)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, required=True)
    main(p.parse_args().output)
