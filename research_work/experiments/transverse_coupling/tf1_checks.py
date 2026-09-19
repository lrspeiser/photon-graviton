"""Suite job for TF-1: reruns the fast exact gates and anchors them to evidence/exact-v1/results.json.
Exit 0 only if every gate passes and the anchors match."""
import json
import sys
import time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import transverse as T   # noqa: E402
from common import read  # noqa: E402


def main():
    t0 = time.time()
    field = T.load_field()
    axis, phi = field['axis'], field['phi']
    arch = read(HERE/'evidence/exact-v1/results.json')
    steer = T.steering_rays_2d(axis, phi)
    index = T.index_rays_2d(axis, phi)
    cmp = T.compare_rays(steer, index)
    probes = T.probes_on_field(axis, phi)
    g6 = T.lagrangian_family(1.)
    g7 = T.velocity_averages(N=200_000)
    g8 = T.focusing_geometry()
    orbits = T.circular_and_radial()
    gates = dict(field_hash=field['sha256'] == arch['field']['sha256'],
                 G1=T.index_speed_control(axis, phi) > 1e-3 and g8['max_speed_drift'] < 1e-9,
                 G2=probes['steering_difference'] < 1e-14 and abs(probes['cwc1_ratio'] - 2) < 1e-9,
                 G5=cmp['all_inward'] and cmp['relative_difference'] < .02,
                 G6=g6['euler_lagrange_residual'] < 1e-6 and g6['speed_change_fraction'] < 1e-8,
                 G7=g7['support_max_difference'] < 5e-3 and g7['isotropic_exact_error'] < 1e-10,
                 G8=g8['focusing_fraction']['transverse'] > .999 and g8['focusing_fraction']['vector_max_abs'] < .1,
                 G9=orbits['radial_steering']['direction_change'] < 1e-12 and orbits['steering']['radius_drift'] < 1e-8)
    anchors = dict(steering_angles=float(np.max(np.abs(np.array([r['angle'] for r in steer]) - np.array([r['angle'] for r in arch['G5']['rays']])))),
                   focusing_transverse=abs(g8['focusing_fraction']['transverse'] - arch['G8']['focusing_fraction']['transverse']),
                   focusing_vector=abs(g8['focusing_fraction']['vector'] - arch['G8']['focusing_fraction']['vector']),
                   support_factor=float(np.max(np.abs(np.array(g7['support_quadrature']) - np.array(arch['G7']['support_quadrature'])))))
    ok = all(gates.values()) and all(v < 1e-9 for v in anchors.values())
    out = dict(gates=gates, anchors=anchors, anchor_matches_archive=all(v < 1e-9 for v in anchors.values()), passed=ok, seconds=time.time() - t0)
    print(json.dumps(out, indent=1))
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
