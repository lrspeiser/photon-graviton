"""Post-evaluation frame diagnostic, not a corrected catalog or a new fit."""
from pathlib import Path
import json
import hashlib
import numpy as np
from scipy.interpolate import PchipInterpolator

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'milky-way-capture/inputs.json'
NEW=HERE.parent/'cepheid-rotation-check/observations.json'
old=json.loads(OLD.read_text())['eilers']['rows']
new=json.loads(NEW.read_text())['rows']
r=np.array([row['R_kpc'] for row in old]);v=np.array([row['vc_kms'] for row in old])
rc=np.array([row['R_kpc'] for row in new]);vc=np.array([row['vc_kms'] for row in new])
delta_r=8.275-8.122
delta_v=250.2-245.8
assert np.all(np.diff(r)>0) and min(rc-delta_r)>=min(r) and max(rc)<=max(r)
interp=PchipInterpolator(r,v,extrapolate=False)
before=vc-interp(rc)
after=(vc-delta_v)-interp(rc-delta_r)
linear_before=vc-np.interp(rc,r,v)
linear_after=vc-delta_v-np.interp(rc-delta_r,r,v)
def metrics(x):
    return dict(rms_difference_kms=float(np.sqrt(np.mean(x*x))),
                mean_difference_kms=float(x.mean()),median_difference_kms=float(np.median(x)),
                positive_rows=int(np.sum(x>0)))
rows=[]
for i,row in enumerate(new):
    rows.append(dict(R_Cepheid_published_kpc=float(rc[i]),Cepheid_published_kms=float(vc[i]),
        Eilers_interpolated_same_published_radius_kms=float(interp(rc[i])),
        published_curve_difference_kms=float(before[i]),
        central_ray_radius_in_old_frame_kpc=float(rc[i]-delta_r),
        illustrative_Cepheid_speed_in_old_frame_kms=float(vc[i]-delta_v),
        central_ray_difference_kms=float(after[i])))
# Check the coordinate identity on the Sun-Galactic-center line at fixed
# heliocentric x. This does not reconstruct a Jeans circular-speed estimator.
heliocentric_x=rc-8.275
reconstructed_r=8.122+heliocentric_x
assert np.max(abs(reconstructed_r-(rc-delta_r)))<1e-12
result=dict(classification='Post-exposure central-ray frame sensitivity; not a revised observed curve, exact Jeans correction or fresh validation',
    published_solar_frames=dict(Eilers=dict(R_kpc=8.122,V_phi_kms=245.8),Feng=dict(R_kpc=8.275,V_phi_kms=250.2)),
    delta_R_kpc=delta_r,delta_V_phi_kms=delta_v,
    illustrative_transform='R_old = R_Cepheid - delta_R; vphi_old = vphi_Cepheid - delta_V_phi, on central ray with unchanged heliocentric position/velocity. Applied to Vc only as a constant-drift approximation.',
    unadjusted_curve_difference=metrics(before),central_ray_constant_drift_difference=metrics(after),
    interpolation_sensitivity=dict(max_PCHIP_linear_before_kms=float(max(abs(before-linear_before))),
        max_PCHIP_linear_after_kms=float(max(abs(after-linear_after)))),
    velocity_projection_only_at_fixed_axes=dict(azimuth_wedge_deg=30,
        delta_vphi_min_kms=float(delta_v*np.cos(np.pi/6)),delta_vphi_max_kms=delta_v),
    published_data_or_previous_scores_modified=False,
    limitations=['Vphi and circular speed Vc are different quantities.',
                 'Changed radial moments, gradients, bin membership, distances and direction bases require individual stellar data.',
                 'Interpolation sensitivity is not observational covariance or a statistical significance.',
                 'The fixed-axis velocity projection range is not a bound on total Vc-estimator systematics.'],
    input_sha256={str(p.relative_to(HERE.parents[2])):hashlib.sha256(p.read_bytes()).hexdigest() for p in [OLD,NEW]})
for name,obj in [('results.json',result),('curve-differences.json',rows)]:
    (HERE/name).write_text(json.dumps(obj,indent=2)+'\n',newline='\n')
print(json.dumps(result,indent=2))
