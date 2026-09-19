"""Read actual archived measurements, then evaluate required prediction capabilities."""
import io,zipfile
import numpy as np
from common import HERE,ROOT,hashes,read

def run(out):
    paths=[ROOT/'temporal_candidate_audit/data/spectral_aging.json',
           ROOT/'temporal_candidate_audit/data/Rotmod_LTG.zip',
           ROOT/'research_work/results/cluster-observation-readiness/kubo-figure-data.json',
           ROOT/'research_work/results/slacs-resolved-input-audit/results.json',
           ROOT/'research_work/results/path-memory/cl2-inputs-xcop-profiles.json']
    source_hashes=hashes(paths)
    aging=read(paths[0]);arr=np.array([r[1:] for r in aging],float)
    data_checks=dict(timing_rows_finite=bool(np.all(np.isfinite(arr))),
                     timing_errors_positive=bool(np.all(arr[:,2]>0)))
    disks=[];points=0
    with zipfile.ZipFile(paths[1]) as z:
        for name in sorted(z.namelist()):
            if not name.endswith('_rotmod.dat'):continue
            values=np.loadtxt(io.BytesIO(z.read(name)))
            if not np.all(np.isfinite(values)):raise ValueError(name)
            disks.append(name);points+=len(values)
    coma=read(paths[2]);cs=coma['rows']
    data_checks['coma_errors_positive']=all(r['plotted_sigma_t']>0 and r['plotted_sigma_cross']>0 for r in cs)
    wanted=['J0037-0942','J1112+0826','J1204+0358','J1402+6321','J1621+3931','J1630+4520']
    systems={r['Name']:r for r in read(paths[3])['systems']}
    lens_rows=[]
    for name in wanted:
        row=systems[name];cov=np.array(row['covariance_kms_squared'])
        eig=float(np.linalg.eigvalsh(cov).min())
        lens_rows.append(dict(name=name,bins=row['n_bins'],minimum_covariance_eigenvalue=eig))
    data_checks['lens_covariance_positive']=all(r['minimum_covariance_eigenvalue']>0 for r in lens_rows)
    # Only ordinary density shells are consumed. No hydro_mass/hydro_params/NFW fit is used.
    clusters=[]
    for name,row in read(paths[4])['clusters'].items():
        d=row['density'];inner=np.array(d['r_in_kpc']);outer=np.array(d['r_out_kpc'])
        ne=np.array(d['ne_cm3'])
        clusters.append(dict(name=name,density_shells=len(ne),valid=bool(np.all(outer>inner)&np.all(ne>=0))))
    data_checks['cluster_density_shells_valid']=all(r['valid'] for r in clusters)
    prerequisites=[
        dict(name='microscopic_identity',ready=False,evidence='Scalar is effective, not a derived particle/spin or charge interaction.'),
        dict(name='matter_rods_clocks',ready=False,evidence='Clock actions are postulated; rods are fixed and no microscopic bound material is solved.'),
        dict(name='three_dimensional_gravity',ready=False,evidence='Only homogeneous, 1D and 2D TE fixed-space Hamiltonians; no universal metric/stress equation.'),
        dict(name='dimensional_normalization',ready=False,evidence='Dimensionless fixtures have no derived universal conversion to stellar or cluster units.'),
        dict(name='bolometric_source_history',ready=False,evidence='Prepared finite pulses, no inferred light history or formation calculation for actual objects.'),
        dict(name='shared_motion_lensing',ready=False,evidence='Material force depends on I/M; rest-coupling alternatives add a new source and remain postulates.')
    ]
    inventory=[
        dict(target='spectral aging and event timing',count=len(aging),
             available='Published exposed spectral-aging rows; not a fresh holdout.',
             candidate='Toy measured frequency/arrival relation exists.',
             missing='Object-specific source/path/clock evolution and line formation.'),
        dict(target='disk rotation',count=len(disks),points=points,
             available='SPARC published rotation/ordinary-light tables.',
             candidate='Toy material force only.',missing='3D source history, dimensional scale, stellar dynamics.'),
        dict(target='resolved lens motions',count=len(lens_rows),systems=lens_rows,
             available='Archived exposed six-lens kinematic bins and covariances.',
             candidate='No aperture or orbit prediction.',missing='Universal matter completion and fixed-distance 3D prediction.'),
        dict(target='Coma shear',count=len(cs),available=coma['status'],
             candidate='Frozen 2D ray bend only.',missing='3D generated profile, source geometry and observational covariance.'),
        dict(target='cluster ordinary gas',count=len(clusters),clusters=clusters,
             available='X-COP derived ordinary electron-density shells at adopted published radii.',
             candidate='No cluster hydrostatic prediction.',missing='3D source and pressure evolution; physical calibration.')
    ]
    if hashes(paths)!=source_hashes:raise RuntimeError('Input changed during audit')
    return dict(stage='readiness',input_sha256=source_hashes,inventory=inventory,gates=data_checks,
                numerical_pass=all(data_checks.values()),prerequisites=prerequisites,
                observational_ready=all(r['ready'] for r in prerequisites),
                astrophysical_fits_executed=False,
                reason='Required predictions do not exist; no fictional observational chi-square is computed.',
                scope='Static tests only. Published units/radii retained without refitting distances; no halo fields evaluated.')
