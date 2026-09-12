"""Age/size-free conditional deductions; no observational refitting.

Run from any directory. Uses known linear transport/ODE mathematics.
These checks verify deductions from postulates, not a microscopic theory.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import numpy as np
from scipy.integrate import quad, solve_ivp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT / "companion_causal_test/energy.csv"
checks = []


def check(name, error, tolerance):
    checks.append(dict(name=name, error=float(error), tolerance=tolerance,
                       passed=bool(error <= tolerance)))


# Free illumination duration: infer conditional requirements, not cosmic age.
with SOURCE.open(encoding="utf-8", newline="") as f:
    rows = list(csv.DictReader(f))
requirements = []
for row in rows:
    ratio = float(row["full_conversion_shortfall"])
    requirements.append(dict(galaxy=row["galaxy"],
                             required_full_capture_years=ratio * 1e10,
                             equivalent_radiated_solar_masses=float(row["extra_mass_Msun"])))
times = np.array([r["required_full_capture_years"] for r in requirements])
# Check archived ratio against independent physical units; archive constants
# may differ slightly. This is a consistency check, not a new mass inference.
for row, req in zip(rows, requirements):
    direct = float(row["extra_mass_Msun"]) * 1.98847e30 * 299792458.0**2
    direct /= float(row["Lbol_Lsun"]) * 3.828e26 * (365.25 * 86400)
    check("archived energy units " + row["galaxy"],
          abs(direct / req["required_full_capture_years"] - 1), 0.002)

# Time-only changes rescale a fixed spatial source, not its normalized shape.
profile = np.array([0.2, 0.7, 1.3, 0.9])
for duration in [1e-8, 1, 1e8]:
    scaled = duration * profile
    check("normalized shape duration " + str(duration),
          np.max(abs(scaled / scaled.sum() - profile / profile.sum())), 1e-14)

# Dimensionless homogeneous energy exchange: gamma -> companion -> deposit.
# Parameters are diagnostic values, NOT chosen cosmic times or inferred rates.
reservoir_cases = []
for a, b in [(0.3, 0.7), (1, 1), (2, 0.2)]:
    for finite_fuel in [False, True]:
        fuel_time = 2.0
        def rhs(t, y):
            j = math.exp(-t / fuel_time) if finite_fuel else 1.0
            return [j-a*y[0], a*y[0]-b*y[1], b*y[1]]
        end = 200 / min(a, b)
        ts = np.linspace(0, end, 501)
        sol = solve_ivp(rhs, [0, end], [0, 0, 0], t_eval=ts,
                        rtol=2e-11, atol=2e-12)
        assert sol.success
        total = fuel_time * (-np.expm1(-ts / fuel_time)) if finite_fuel else ts
        error = np.max(abs(sol.y.sum(axis=0)-total)) / max(1, total[-1])
        label = f"a={a} b={b} finite_fuel={finite_fuel}"
        check("energy conservation " + label, error, 1e-9)
        check("nonnegative energy " + label, max(0, -sol.y.min()), 1e-9)
        if finite_fuel:
            check("finite fuel deposit limit " + label,
                  abs(sol.y[2, -1] / fuel_time - 1), 1e-8)
        else:
            check("permanent deposit late growth " + label,
                  abs(b*sol.y[1, -1]-1), 1e-8)
        reservoir_cases.append(dict(a=a, b=b, finite_fuel=finite_fuel,
                                    final_energy=sol.y[:, -1].tolist(),
                                    final_time=end))

# Constant deposited power with finite lifetime: U=P*tau*(1-exp(-T/tau)).
# Lost deposit energy is put into a fourth, explicit receiver, never erased.
for tau in [0.2, 1, 20]:
    ts = np.linspace(0, 20*tau, 101)
    sol = solve_ivp(lambda t,y: [1-y[0]/tau, y[0]/tau],
                    [0, ts[-1]], [0,0], t_eval=ts, rtol=1e-11, atol=1e-12)
    assert sol.success
    analytic = tau * (-np.expm1(-ts/tau))
    check("finite lifetime solution " + str(tau),
          np.max(abs(sol.y[0]-analytic))/tau, 1e-9)
    check("finite lifetime receiver ledger " + str(tau),
          np.max(abs(sol.y.sum(axis=0)-ts))/ts[-1], 1e-9)

# Infinite-radius shell integrals: constant Euclidean emissivity, steady
# sources, common constant speed, no event stretching, linear capture opacity.
# The capture opacity is a coarse-grained well-interception toy, not void storage.
shells = []
for a, b in [(1, 0), (1, 1), (1, 0.2), (1, 3)]:
    def photon(r):
        return math.exp(-a*r)
    def companion(r):
        if b == 0:
            return -math.expm1(-a*r)
        if a == b:
            return a*r*math.exp(-a*r)
        return a/(b-a)*(math.exp(-a*r)-math.exp(-b*r))
    for radius in [0.1, 1, 10, 100]:
        ug = quad(photon, 0, radius, epsabs=1e-11)[0]
        uc = quad(companion, 0, radius, epsabs=1e-11)[0]
        if b == 0:
            check(f"lossless ray energy R={radius}", abs(ug+uc-radius), 1e-9)
        for r in np.linspace(0, radius, 41):
            pg, pc = photon(r), companion(r)
            check(f"ray positivity {a},{b},{radius},{r}",
                  max(0, -pg, -pc, pg+pc-1), 1e-12)
        shells.append(dict(alpha=a, beta=b, radius=radius,
                           photon_integral=ug, companion_integral=uc))
    if b > 0:
        ug = quad(photon, 0, np.inf, epsabs=1e-11)[0]
        uc = quad(companion, 0, np.inf, epsabs=1e-11)[0]
        check(f"infinite photon integral {b}", abs(ug-1/a), 1e-9)
        check(f"infinite companion integral {b}", abs(uc-1/b), 1e-9)

# Stationary travel law: even enormous constant delays preserve event spacing.
# Analytic derivative is exactly 1; these are well-resolved arithmetic probes.
for delay in [0, 1, 1e6]:
    separation = ((delay+12)-delay)
    check("stationary event separation " + str(delay), abs(separation-12), 0)

out = dict(
    scope="Conditional age/size-free transport and energy deductions; no new observational fit",
    source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    age_requirement_years=dict(median=float(np.median(times)),
                              p10=float(np.quantile(times,.1)),
                              p90=float(np.quantile(times,.9))),
    per_galaxy_requirements=requirements,
    reservoir_cases=reservoir_cases, shell_cases=shells,
    checks=checks, all_checks_passed=all(c["passed"] for c in checks))
(HERE / "results.json").write_text(json.dumps(out, indent=2)+"\n", encoding="utf-8", newline="\n")
print(json.dumps({"checks":len(checks), "passed":out["all_checks_passed"],
                  "age_requirement_years":out["age_requirement_years"],
                  "failures":[c for c in checks if not c["passed"]]}, indent=2))
if not out["all_checks_passed"]:
    raise SystemExit(1)
