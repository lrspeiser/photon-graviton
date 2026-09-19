"""E1: angular radiation activation using ordinary-light morphology proxies."""
import io
import zipfile
import numpy as np
from scipy.special import ndtri
from scipy.stats import qmc
from common import ROOT, PM, read_json, track, inv_cdf, directions

RADII = np.array([0., .5, 1., 2.])

def activation_tensor(P):
    u = np.trace(P)
    return 0. if u == 0 else 27*np.linalg.det(P)/u**3

def activation(points, scale, soft):
    result = []
    for radius in RADII:
        delta = np.array([radius*scale, 0., 0.])-points
        d2 = np.sum(delta*delta, axis=1)
        n = np.divide(delta, np.sqrt(d2)[:, None], out=np.zeros_like(delta),
                      where=d2[:, None] > 0)
        P = (n.T/(d2+(soft*scale)**2)) @ n
        result.append(activation_tensor(P))
    return np.array(result)

def controls():
    rng = np.random.default_rng(1909)
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    P = np.diag([.2, .3, .5])
    expected = np.array([1., 0., 0., 0.])
    got = np.array([activation_tensor(np.eye(3)), activation_tensor(np.diag([1.,0,0])),
                    activation_tensor(np.diag([.5,.5,0])), activation_tensor(np.zeros((3,3)))])
    return dict(exact_error=np.max(abs(got-expected)),
                rotation_error=abs(activation_tensor(P)-activation_tensor(Q@P@Q.T)),
                subdivision_error=abs(activation_tensor(P)-activation_tensor(sum([P/17]*17))),
                luminosity_error=max(abs(activation_tensor(P)-activation_tensor(P*f))
                                     for f in [.001, 1000.]))

def run():
    seq = qmc.Sobol(4, scramble=True, seed=1909).random_base2(14)
    catalogue = {}
    for line in track(ROOT/'temporal_candidate_audit/data/SPARC_Lelli2016c.mrt').read_text().splitlines():
        f = line.split()
        if len(f) == 19:
            try: catalogue[f[0]] = float(f[11])
            except ValueError: pass
    disks, omitted = [], []
    with zipfile.ZipFile(track(ROOT/'temporal_candidate_audit/data/Rotmod_LTG.zip')) as z:
        for name in sorted(z.namelist()):
            if not name.endswith('_rotmod.dat'): continue
            key = name.rsplit('_rotmod', 1)[0]
            data = np.loadtxt(io.BytesIO(z.read(name)))
            rd = catalogue.get(key, 0.)
            good = data[:,6] > 0
            if rd <= 0 or np.count_nonzero(good) < 2:
                omitted.append(key)
                continue
            rr, sb = data[good,0], data[good,6]
            r = np.linspace(0, rr[-1], 2048)
            lum = np.interp(r, rr, sb)
            sample_r = inv_cdf(r, r*lum, seq[:,0])
            phi = 2*np.pi*seq[:,1]
            for h in [.05, .1, .2]:
                points = np.column_stack((sample_r*np.cos(phi), sample_r*np.sin(phi),
                                          h*rd*ndtri(seq[:,2])))
                row = dict(name=key, height=h, scale_kpc=rd,
                           scores={str(n): activation(points[:n], rd, .05)
                                   for n in [4096,8192,16384]},
                           softening_scores={str(s): activation(points,rd,s) for s in [.025,.1]})
                disks.append(row)
    clusters = []
    for name, c in read_json(PM/'cl2-inputs-xcop-profiles.json')['clusters'].items():
        scale = c['header']['R500_kpc']
        d = c['density']
        mid = (np.array(d['r_in_kpc'])+np.array(d['r_out_kpc']))/2
        r = np.linspace(0, d['r_out_kpc'][-1], 2048)
        ne = np.interp(r, mid, d['ne_cm3'])
        p = c['pressure']['xray']
        pressure = np.interp(r, np.array(p['r_over_R500'])*scale, p['P_over_P500'])
        j = ne**2*np.sqrt(pressure/ne)
        radius = inv_cdf(r,r*r*j,seq[:,0])
        points = radius[:,None]*directions(seq[:,1],seq[:,2])
        clusters.append(dict(name=name,scale_kpc=scale,
                             scores={str(n):activation(points[:n],scale,.05)
                                     for n in [4096,8192,16384]},
                             softening_scores={str(s):activation(points,scale,s) for s in [.025,.1]}))
    comparisons=[]
    for soft in ['0.025','0.05','0.1']:
        def score(row):
            return row['scores']['16384'] if soft=='0.05' else row['softening_scores'][soft]
        cc = np.array([score(c) for c in clusters])
        for h in [.05,.1,.2]:
            dd = np.array([score(d) for d in disks if d['height']==h])
            comparisons.append(dict(softening=soft,height=h,
                                    cluster_median=np.median(cc,axis=0),
                                    disk_median=np.median(dd,axis=0),
                                    ratio=np.median(cc,axis=0)/np.median(dd,axis=0),
                                    cluster_10_90=np.quantile(cc,[.1,.9],axis=0),
                                    disk_10_90=np.quantile(dd,[.1,.9],axis=0)))
    allrows=disks+clusters
    change=np.median([np.max(abs(x['scores']['16384']-x['scores']['8192'])) for x in allrows])
    vals=np.concatenate([x['scores'][str(n)] for x in allrows for n in [4096,8192,16384]])
    ctl=controls()
    gates=dict(controls=max(ctl.values())<1e-10,bounds=vals.min()>=-1e-10 and vals.max()<=1+1e-10,
               convergence=change<.05)
    return dict(experiment='E1',controls=ctl,gates=gates,numerical_pass=all(gates.values()),
                radii_over_scale=RADII,median_resolution_change=change,
                disk_count=len(disks)//3,cluster_count=len(clusters),omitted=omitted,
                mechanism_pass=all(np.all(x['ratio']>=3) for x in comparisons),
                comparisons=comparisons,disks=disks,clusters=clusters,
                observational_status='Unresolved: shape proxies, disk-only light, finite softening; no bolometric maps.')
