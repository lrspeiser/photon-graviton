"""Run GF-1 declared screens. Output directories are immutable: no overwrite."""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time
import numpy as np
from model import Mechanics, catalogue, cross, initial, leader_angle, rk4

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]

def clean(value):
    if isinstance(value, np.ndarray): return clean(value.tolist())
    if isinstance(value, (np.integer,)): return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, dict): return {k: clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [clean(v) for v in value]
    if isinstance(value, np.bool_): return bool(value)
    return value

def write_json(path, value):
    path.write_text(json.dumps(clean(value), indent=2, allow_nan=False)+"\n", encoding="utf8")

def manifest():
    return dict(python=platform.python_version(), numpy=np.__version__,
                platform=platform.platform(), command=" ".join(__import__("sys").argv),
                git_head=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
                source_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in sorted(HERE.iterdir()) if p.suffix in (".py", ".md", ".json")},
                prohibited_components={"dark_matter": False, "expanding_background": False,
                                       "fitted_distances": False})

def diagnostics(mech, state, t, chain, amplitude):
    z = mech.terms(state)
    x, v = state[..., :2], state[..., 2:4]
    xc = x-np.mean(x, axis=1, keepdims=True)
    vc = v-np.mean(v, axis=1, keepdims=True)
    radii = np.linalg.norm(xc, axis=-1)
    denom = radii*np.linalg.norm(vc, axis=-1)
    circ_each = np.divide(cross(xc, vc), denom, out=np.zeros_like(denom), where=denom != 0)
    rr = z["r"].copy()
    rr[:, np.arange(x.shape[1]), np.arange(x.shape[1])] = np.inf
    reaction = mech.rhs(t, state, chain, amplitude)[1] if chain else np.zeros((len(state), 2))
    return dict(energy=z["energy"], energy_scale=.5*np.sum(v*z["p"], axis=(1, 2))+np.abs(z["U"]),
                momentum=np.sum(z["p"], axis=1),
                angular=np.sum(cross(x, z["p"]), axis=1),
                angular_scale=np.sum(np.abs(cross(x, z["p"])), axis=1),
                speed_scale=np.sum(np.linalg.norm(v, axis=-1), axis=1),
                radius=np.mean(radii, axis=1),
                spread=np.std(radii, axis=1)/np.mean(radii, axis=1),
                circulation=np.abs(np.mean(circ_each, axis=1)),
                separation=np.min(rr, axis=(1, 2)),
                guide_work=z["guide_work"],
                guide_force=np.sum(z["guide"], axis=1),
                reaction=reaction, external_power=np.sum(reaction*v[:, 0], axis=1),
                leader_speed=np.linalg.norm(v[:, 0], axis=-1))

def best_correlation(signal, leader):
    best, lag = -1., 0
    # All nonnegative lags through 20 time units; sample cadence 0.2.
    for offset in range(min(100, len(signal)//2)+1):
        a = signal[offset:]
        b = leader[:len(a)]
        if np.std(a) <= 1e-14 or np.std(b) <= 1e-14: continue
        score = np.corrcoef(a, b)[0, 1]
        if score > best: best, lag = float(score), offset*.2
    return best, lag

def summarize(row, history, path, period, fixture, amplitude, duration, first_failure, min_step_sep):
    h = history
    v = path[:, :, 2:4]
    result = dict(id=row["id"], fixture=fixture, duration=duration,
                  initial_period=period, periods_covered=duration/period,
                  first_nonfinite_time=first_failure,
                  min_separation=min_step_sep,
                  circulation_min=float(np.min(h["circulation"])),
                  circulation_final=float(h["circulation"][-1]),
                  radius_ratio_min=float(np.min(h["radius"])/h["radius"][0]),
                  radius_ratio_max=float(np.max(h["radius"])/h["radius"][0]),
                  radius_ratio_final=float(h["radius"][-1]/h["radius"][0]),
                  relative_radius_spread_max=float(np.max(h["spread"])),
                  energy_relative_drift=float(np.max(np.abs(h["energy"]-h["energy"][0]))/h["energy_scale"][0]),
                  momentum_relative_drift=float(np.max(np.linalg.norm(h["momentum"]-h["momentum"][0], axis=-1))/max(h["speed_scale"][0], 1e-30)),
                  angular_relative_drift=float(np.max(np.abs(h["angular"]-h["angular"][0]))/max(h["angular_scale"][0], 1e-30)),
                  guide_work_max=float(np.max(np.abs(h["guide_work"]))),
                  guide_force_max=float(np.max(np.linalg.norm(h["guide_force"], axis=-1))),
                  external_work_sampled=float(np.trapezoid(h["external_power"], h["time"])),
                  leader_speed_error=float(np.max(np.abs(h["leader_speed"]-1))) if fixture=="chain" else None)
    if fixture == "chain":
        angle = np.arctan2(v[..., 1], v[..., 0])
        gains = np.max(np.abs(angle[:, 1:]), axis=0)/max(abs(amplitude), 1e-30)
        sig = np.mean(angle[:, -4:], axis=1)
        lead = np.array([leader_angle(t, amplitude) for t in h["time"]])
        corr, lag = best_correlation(sig, lead)
        arrival = [float(h["time"][np.flatnonzero(np.abs(angle[:, i]) > .02)[0]])
                   if np.any(np.abs(angle[:, i]) > .02) else None for i in range(1, angle.shape[1])]
        result.update(chain_gain=float(np.mean(gains[-4:])), maximum_follower_gain=float(np.max(gains)),
                      gain_per_follower=gains.tolist(), arrival_time_per_follower=arrival,
                      last_four_correlation=corr, best_lag=lag,
                      no_turn_transverse_max=float(np.max(np.abs(v[..., 1]))) if amplitude == 0 else None)
        result["behavioral_pass"] = bool(amplitude != 0 and result["chain_gain"] >= .1
             and corr >= .7 and result["maximum_follower_gain"] <= 3 and min_step_sep > .1
             and first_failure is None)
    else:
        result["behavioral_pass"] = bool(result["circulation_min"] >= .8
             and result["radius_ratio_min"] >= .7 and result["radius_ratio_max"] <= 1.3
             and result["relative_radius_spread_max"] < .25 and min_step_sep > .1
             and first_failure is None)
    result["conservation_pass"] = bool(row["family"] == 5 and fixture == "ring"
           and result["energy_relative_drift"] < 1e-3 and result["momentum_relative_drift"] < 1e-5
           and result["angular_relative_drift"] < 1e-3 and first_failure is None)
    return result

def simulate(rows, fixture, out, dt=.04, duration=None, n=None, amplitude=.4, coupling=.2):
    out.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter()
    mech = Mechanics(rows, coupling=coupling)
    y, periods = initial(mech, fixture, n)
    duration = duration or (32 if fixture == "chain" else 64)
    steps = int(round(duration/dt))
    cadence = int(round(.2/dt))
    if abs(steps*dt-duration) > 1e-12 or abs(cadence*dt-.2)>1e-12:
        raise ValueError("dt must divide duration and record cadence")
    samples, paths, times = [], [], []
    failures = [None]*len(rows)
    min_sep = np.full(len(rows), np.inf)
    live = np.ones(len(rows), dtype=bool)
    dummy = y.copy()
    for step in range(steps+1):
        t = step*dt
        # Closest approach on every integration step, not just saved frames.
        d = y[:, None, :, :2]-y[:, :, None, :2]
        sep = np.linalg.norm(d, axis=-1)
        sep[:, np.arange(y.shape[1]), np.arange(y.shape[1])] = np.inf
        min_sep = np.minimum(min_sep, np.min(sep, axis=(1, 2)))
        if step % cadence == 0:
            diag = diagnostics(mech, y, t, fixture == "chain", amplitude)
            for a in diag.values(): a[~live] = np.nan
            samples.append(diag)
            p = y.copy(); p[~live] = np.nan
            paths.append(p)
            times.append(t)
        if step == steps: break
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            y = rk4(mech, t, y, dt, fixture == "chain", amplitude)
        bad = ~np.all(np.isfinite(y), axis=(1, 2))
        for i in np.flatnonzero(bad & live): failures[i] = (step+1)*dt
        live &= ~bad
        # Inactive batch rows are masked out of all evidence and all pass decisions.
        # This prevents a failed row from overflowing unrelated batch members.
        y[~live] = dummy[~live]
    hist = {key: np.stack([s[key] for s in samples], axis=1) for key in samples[0]}
    paths = np.stack(paths, axis=1)
    hist["time"] = np.array(times)
    summaries = []
    for i, row in enumerate(rows):
        individual = {k: (a if k=="time" else a[i]) for k, a in hist.items()}
        summaries.append(summarize(row, individual, paths[i], periods[i], fixture,
                                   amplitude, duration, failures[i], min_sep[i]))
    np.savez_compressed(out/"trajectories.npz", state=paths,
                        ids=np.array([r["id"] for r in rows]), **hist)
    write_json(out/"summary.json", dict(fixture=fixture, dt=dt, duration=duration, n=y.shape[1],
              amplitude=amplitude, coupling=coupling, wall_seconds=time.perf_counter()-start,
              failures=failures, results=summaries))
    print(f"{out.name}: {len(rows)} {fixture} cases; "
          f"{sum(r['behavioral_pass'] for r in summaries)} behavioral passes; "
          f"{time.perf_counter()-start:.1f}s", flush=True)
    return summaries

def combined_screen(root):
    rows = catalogue()
    results = {}
    for label, batch in [("phenomenological", rows[:600]), ("conservative", rows[600:])]:
        results[label] = {}
        for fixture in ("chain", "ring"):
            results[label][fixture] = simulate(batch, fixture, root/f"{label}-{fixture}")
    write_json(root/"screen.json", results)
    return results

def select(screen, group, count):
    a = {r["id"]: r for r in screen[group]["chain"]}
    b = {r["id"]: r for r in screen[group]["ring"]}
    def rank(identifier):
        gain = a[identifier]["chain_gain"]
        circ = b[identifier]["circulation_min"]
        return (-int(a[identifier]["behavioral_pass"])-int(b[identifier]["behavioral_pass"]),
                -min(1, gain if gain is not None else -1, circ if circ is not None else -1), identifier)
    return sorted(a, key=rank)[:count]

def refinement(screen, root):
    allrows = {r["id"]: r for r in catalogue()}
    selected = {g: select(screen, g, n) for g, n in (("phenomenological", 12), ("conservative", 3))}
    write_json(root/"selection.json", selected)
    runs = {}
    for group, ids in selected.items():
        batch = [allrows[i] for i in ids]
        for fixture in ("chain", "ring"):
            name = f"{group}-{fixture}-dt002"
            runs[name] = simulate(batch, fixture, root/name, dt=.02)
        if group == "phenomenological":
            for name, args in [("small-turn", dict(amplitude=.2)),
                               ("large-turn", dict(amplitude=.6)),
                               ("reverse-turn", dict(amplitude=-.4)),
                               ("long-chain", dict(n=24)),
                               ("reverse-coupling", dict(coupling=-.2)),
                               ("zero-input", dict(amplitude=0.))]:
                runs[name] = simulate(batch, "chain", root/name, **args)
            runs["reverse-coupling-ring"] = simulate(batch, "ring", root/"reverse-coupling-ring", coupling=-.2)
        else:
            runs["conservative-long-ring"] = simulate(batch, "ring", root/"conservative-long-ring", dt=.02, duration=256)
            runs["conservative-zero-input"] = simulate(batch, "chain", root/"conservative-zero-input", amplitude=0.)
    write_json(root/"refinement.json", runs)
    return runs

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("phase", choices=["catalogue", "screen", "refine"])
    parser.add_argument("--output", type=Path)
    parser.add_argument("--screen", type=Path)
    args = parser.parse_args()
    if args.phase == "catalogue":
        rows = catalogue()
        write_json(HERE/"catalogue.json", rows)
        text = ["# GF-1: 618 specified candidate equations\n",
                "These are variants and parameter cases of credited mathematical constructions, not 618 novelty claims.",
                "Definitions, units, shared U, initial conditions and zero-speed conventions: [protocol](protocol.md).\n"]
        for r in rows:
            text += [f"## {r['id']}\n", chr(96)+r["equation"]+chr(96)+"\n", r["status"]+"\n"]
        # Use literal Markdown code spans without escaping backticks.
        (HERE/"catalogue.md").write_text("\n".join(text), encoding="utf8")
        print(f"Wrote {len(rows)} unique formulas", flush=True)
        return
    if args.output is None: parser.error("--output required")
    args.output.mkdir(parents=True, exist_ok=False)
    write_json(args.output/"manifest.json", manifest())
    if args.phase == "screen": combined_screen(args.output)
    else:
        if args.screen is None: parser.error("--screen required")
        refinement(json.loads(args.screen.read_text(encoding="utf8")), args.output)
    write_json(args.output/"complete.json", dict(complete=True))

if __name__ == "__main__":
    main()
