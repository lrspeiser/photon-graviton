"""Compare proposed force signals with one measured truncation sensitivity."""
import json
import numpy as np
from run import FOUNDATION,FIELD_CACHE,save
from fast_multipole import FastMultipole
from cached_companion import CachedCompanion

base=[r for r in json.loads((FOUNDATION/'field-predictions.json').read_text()) if r['model']=='ordinary_matter']
xyz=np.array([[r['x_kpc'],r['y_kpc'],r['z_kpc']] for r in base])
_,low=FastMultipole.load(FIELD_CACHE/'bar-L40.npz').evaluate(xyz)
_,high=FastMultipole.load(FIELD_CACHE/'bar-L64.npz').evaluate(xyz)
shift=np.linalg.norm(high-low,axis=1)
rows=[]
for geometry in ['equatorial','caps','shell']:
    _,force=CachedCompanion(geometry).evaluate(xyz)
    signal=np.linalg.norm(force,axis=1)
    ratio=signal/np.maximum(shift,1e-30)
    rows.append(dict(geometry=geometry,positions=72,
        signal_below_bar_resolution_change=int(np.sum(signal<shift)),
        minimum_signal_to_bar_resolution_change=float(ratio.min()),
        median_signal_to_bar_resolution_change=float(np.median(ratio)),
        minimum_ratio_position_kpc=xyz[ratio.argmin()].tolist()))
save('effect-resolution.json',dict(comparisons=rows,
    interpretation='A finite bar-order difference is a sensitivity indicator, not a rigorous force-error bound. Does not include ordinary-matter parameter uncertainty, disk/nuclear convergence or observational errors.',
    observed_data_tested=False))
print(rows)
