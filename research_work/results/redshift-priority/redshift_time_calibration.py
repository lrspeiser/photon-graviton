"""Observed stretch targets, not local-clock measurements or fitted physics."""
import csv
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
source = ROOT / 'redshift_paper/all_164_groups.csv'
rows = list(csv.DictReader(source.open(encoding='utf-8-sig')))
fields = ['pgc', 'group_pgc', 'distance_mpc', 'observed_cmb_z',
          'observed_stretch_factor', 'log_observed_stretch',
          'all_shift_time_average_rate_per_mpc', 'all_shift_time_photon_energy_fraction',
          'effective_signal_clock_ratio', 'void_rate_if_fraction_1_per_mpc',
          'void_rate_if_fraction_0p9_per_mpc', 'void_rate_if_fraction_0p75_per_mpc']
with (HERE / 'observed-time-stretch-targets.csv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, lineterminator='\n')
    writer.writeheader()
    for row in rows:
        z = float(row['observed_cmb_z'])
        distance = float(row['catalog_distance_mpc'])
        assert z > -1 and distance > 0
        stretch = 1 + z
        depth = math.log1p(z)
        assert math.isclose(math.expm1(depth), z, rel_tol=1e-13, abs_tol=1e-15)
        writer.writerow(dict(zip(fields, [row['pgc'], row['group_pgc'], distance, z,
            stretch, depth, depth/distance, z/stretch, 1/stretch,
            depth/distance, depth/(0.9*distance), depth/(0.75*distance)])))
assert len(rows) == 164 and len({r['group_pgc'] for r in rows}) == 164
metadata = {
    'source': source.relative_to(ROOT).as_posix(),
    'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
    'rows': len(rows),
    'status': 'Previously exposed observed-redshift transformations; no fit or fresh validation.',
    'assumptions': ['Published distances retained as stipulated facts.',
                    'Observed CMB z kept unchanged; no nuisance correction invented.',
                    'Columns prefixed all_shift_time assume the entire observed shift is due to the time mechanism.',
                    'Per-object averages are diagnostic targets, never independently adjustable prediction parameters.',
                    'Energy fraction additionally assumes E proportional to frequency and all modeled loss enters companions.',
                    'None of these columns measures a physical local clock rate.'],
    'formula_provenance': 'Established redshift definition, logarithm and photon energy relation; conditional algebraic transformations here, not novel laws.'
}
(HERE / 'observed-time-stretch-targets.json').write_text(json.dumps(metadata, indent=2) + '\n', newline='\n')
print(json.dumps({'rows': len(rows), 'source_sha256': metadata['source_sha256']}))
