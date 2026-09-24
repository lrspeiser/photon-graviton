"""Round 16, step 2a: what opens the quiet store? The velocity of matter relative to the companion flowing past it.

The independent calculation (independent-r15/) mixes a piece's quiet mode D into its radiating modes B_m by
delta_m = chi w_m, with w a 'relative velocity', and needs a separate quiet leak gamma_0; its heat weight is then
k = 3 chi^2 sigma^2/(gamma_0 gamma) = 3 sigma^2/u_eff^2 with u_eff^2 = gamma_0 gamma/chi^2, a relation between unknown rates.
Hypothesis tested here: w is the piece's velocity relative to the companion flowing through it, w = v - V_flow. The
companion streams away from the matter that emits it at the speed u, and keeps its sources' mean motion (round 4's
memory). Then, with NO quiet leak of its own (gamma_0 = 0):
  * matter at rest in its own flow sees |w| = u: its quiet leak is the flow's doing, gamma_0 = chi^2 u^2/gamma, so the
    relation u^2 = gamma_0 gamma/chi^2 holds identically, with u the flow's speed;
  * random motion (1D spread sigma) adds 3 sigma^2: released/cold = 1 + 3 sigma^2/u^2, the law's heat weight, factor 3 and
    all, with the same u as the travel speed;
  * collisions (velocities redrawn at rate nu) scramble only the random part: 1 + (3 sigma^2/u^2) gamma/(gamma + nu);
  * a steady drift along the flow at speed v_r gives (u - v_r)^2 + ...: outflow at u switches the leak off, infall adds.
The equations are the calculation's own (integrated exactly as it writes them, RK4, gamma = 1), with delta_m(t) =
chi (v_m(t) - U s_m): U = chi u/gamma the flow term, s = the flow's direction, v the piece's velocity.

    python code/stream_store_v16.py --output run-stream-store-v16/stream_store_v16.json
"""
from __future__ import annotations
import argparse, json, math, time
from pathlib import Path
import numpy as np
from numba import njit


@njit(cache=True)
def rate(x, y0, y1, y2, a0, a1, a2, g0):
    # the calculation's finite store (x = D, y_m = B_m as real quadratures), gamma = 1
    return (-g0 * x - a0 * y0 - a1 * y1 - a2 * y2, -y0 + a0 * x, -y1 + a1 * x, -y2 + a2 * x,
            2 * g0 * x * x + 2 * (y0 * y0 + y1 * y1 + y2 * y2))


@njit(cache=True)
def sim(U, q, nu, drift, seed, cells, T, burn, dt, g0):
    """Mean radiated power per unit stored energy over [burn, T], averaged over cells. delta = v - U s with s = z;
    v = drift z + random part (per component rms q; Ornstein-Uhlenbeck at rate nu, or fixed when nu = 0)."""
    np.random.seed(seed)
    rho = math.exp(-nu * dt); kick = q * math.sqrt(max(0.0, 1 - rho * rho))
    nstep = int(round(T / dt)); nb = int(round(burn / dt))
    acc = 0.0
    for j in range(cells):
        v0 = q * np.random.normal(); v1 = q * np.random.normal(); v2 = q * np.random.normal()
        x = 1.0; y0 = 0.0; y1 = 0.0; y2 = 0.0
        num = 0.0; den = 0.0
        for k in range(nstep):
            if nu > 0:
                v0 = rho * v0 + kick * np.random.normal(); v1 = rho * v1 + kick * np.random.normal(); v2 = rho * v2 + kick * np.random.normal()
            a0 = v0; a1 = v1; a2 = v2 + drift - U
            r1 = rate(x, y0, y1, y2, a0, a1, a2, g0)
            r2 = rate(x + .5 * dt * r1[0], y0 + .5 * dt * r1[1], y1 + .5 * dt * r1[2], y2 + .5 * dt * r1[3], a0, a1, a2, g0)
            r3 = rate(x + .5 * dt * r2[0], y0 + .5 * dt * r2[1], y1 + .5 * dt * r2[2], y2 + .5 * dt * r2[3], a0, a1, a2, g0)
            r4 = rate(x + dt * r3[0], y0 + dt * r3[1], y1 + dt * r3[2], y2 + dt * r3[3], a0, a1, a2, g0)
            x += dt * (r1[0] + 2 * r2[0] + 2 * r3[0] + r4[0]) / 6
            y0 += dt * (r1[1] + 2 * r2[1] + 2 * r3[1] + r4[1]) / 6
            y1 += dt * (r1[2] + 2 * r2[2] + 2 * r3[2] + r4[2]) / 6
            y2 += dt * (r1[3] + 2 * r2[3] + 2 * r3[3] + r4[3]) / 6
            if k >= nb:
                num += r1[4]; den += x * x + y0 * y0 + y1 * y1 + y2 * y2
        acc += num / den
    return acc / cells


def group(U, q, nu=0.0, drift=0.0, seeds=8, cells=64, T=60.0, burn=10.0, dt=0.01, g0=0.0):
    v = np.array([sim(U, q, nu, drift, s, cells, T, burn, dt, g0) for s in range(1, seeds + 1)])
    return float(v.mean()), float(v.std(ddof=1) / np.sqrt(len(v)))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    U = 0.01                      # chi u / gamma: the flow's mixing, small (the adiabatic regime of the calculation)
    out = dict(experiment='round 16: the flow opens the store', U=U, gamma=1.0, gamma_0=0.0, rows={})
    cold, cold_e = group(U, 0.0)
    out['rows']['cold'] = dict(leak=cold, err=cold_e, expected=2 * U ** 2)
    print(f'cold (at rest in its own flow, no quiet leak of its own): leak {cold:.4e} per unit store, expected 2U^2 = {2 * U * U:.4e}', flush=True)
    free = []
    for r in (0.25, 0.5, 1.0, 2.0, 4.0, 8.0):          # sigma / u
        m, e = group(U, r * U)
        free.append(dict(sigma_over_u=r, ratio=m / cold, ratio_err=e / cold, law=1 + 3 * r * r))
        print(f'free, sigma = {r:g} u: released/cold {m / cold:.3f} +- {e / cold:.3f}   (1 + 3 sigma^2/u^2 = {1 + 3 * r * r:.3f})', flush=True)
    out['rows']['free'] = free
    coll = []
    r = 2.0
    for nu in (1.0, 3.0, 10.0, 30.0):
        m, e = group(U, r * U, nu=nu)
        k_meas = m / cold - 1; k_free = 3 * r * r
        coll.append(dict(sigma_over_u=r, nu=nu, heat=k_meas, heat_err=e / cold, expected=k_free / (1 + nu), free=k_free))
        print(f'colliding, sigma = 2u, nu = {nu:g} gamma: heat {k_meas:.3f} +- {e / cold:.3f}  (k gamma/(gamma+nu) = {k_free / (1 + nu):.3f}); cold part kept', flush=True)
    out['rows']['collisional'] = coll
    drift = []
    for d in (-2.0, -1.0, -0.5, 0.5, 1.0, 2.0):     # v_r / u, positive along the flow (outflow)
        m, e = group(U, 0.0, drift=d * U)
        drift.append(dict(vr_over_u=d, ratio=m / cold, expected=(1 - d) ** 2))
        print(f'steady drift along the flow v_r = {d:+g} u: leak/cold {m / cold:.3f}  ((1 - v_r/u)^2 = {(1 - d) ** 2:.3f})', flush=True)
    out['rows']['drift'] = drift
    # rotation: a piece moving at v_rot across the flow; with memory the flow co-moves (relative velocity u only),
    # without memory the flow is fixed in the lab: v_rot adds as heat
    rot = []
    for vr in (0.5, 1.0, 2.0):
        m_mem, _ = group(U, 0.0)                                       # co-moving flow: only u
        m_lab = np.mean([sim(U, 0.0, 0.0, 0.0, s, 64, 60.0, 10.0, 0.01, 0.0) for s in range(1, 2)])
        # moving across the flow at vr u: delta = (vr u, 0, -u) -> |w|^2 = u^2 (1 + vr^2); emulate with a fixed x-velocity
        m_x = np.mean([sim_cross(U, vr * U, s) for s in range(1, 5)])
        rot.append(dict(v_rot_over_u=vr, with_memory=m_mem / cold, without_memory=m_x / cold, expected_without=1 + vr * vr))
        print(f'rotation at {vr:g} u: with memory {m_mem / cold:.3f}, without memory {m_x / cold:.3f} (1 + v^2/u^2 = {1 + vr * vr:.3f})', flush=True)
    out['rows']['rotation'] = rot
    out['seconds'] = time.monotonic() - t0
    args.output.write_text(json.dumps(out, indent=1) + '\n')
    print(f'wrote {args.output} ({out["seconds"]:.0f} s)')


@njit(cache=True)
def sim_cross(U, vx, seed):
    """A piece moving steadily across the flow (x) while the flow runs along z: delta = (vx, 0, -U)."""
    dt = 0.01; T = 60.0; burn = 10.0
    nstep = int(round(T / dt)); nb = int(round(burn / dt))
    x = 1.0; y0 = 0.0; y1 = 0.0; y2 = 0.0; num = 0.0; den = 0.0
    for k in range(nstep):
        r1 = rate(x, y0, y1, y2, vx, 0.0, -U, 0.0)
        r2 = rate(x + .5 * dt * r1[0], y0 + .5 * dt * r1[1], y1 + .5 * dt * r1[2], y2 + .5 * dt * r1[3], vx, 0.0, -U, 0.0)
        r3 = rate(x + .5 * dt * r2[0], y0 + .5 * dt * r2[1], y1 + .5 * dt * r2[2], y2 + .5 * dt * r2[3], vx, 0.0, -U, 0.0)
        r4 = rate(x + dt * r3[0], y0 + dt * r3[1], y1 + dt * r3[2], y2 + dt * r3[3], vx, 0.0, -U, 0.0)
        x += dt * (r1[0] + 2 * r2[0] + 2 * r3[0] + r4[0]) / 6
        y0 += dt * (r1[1] + 2 * r2[1] + 2 * r3[1] + r4[1]) / 6
        y1 += dt * (r1[2] + 2 * r2[2] + 2 * r3[2] + r4[2]) / 6
        y2 += dt * (r1[3] + 2 * r2[3] + 2 * r3[3] + r4[3]) / 6
        if k >= nb:
            num += r1[4]; den += x * x + y0 * y0 + y1 * y1 + y2 * y2
    return num / den


if __name__ == '__main__':
    main()
