"""Describe like-for-like training overlap; no force inference or held-out scores."""
import numpy as np
import pandas as pd
from run import CACHE,H,save,coordinates,region

d=pd.read_parquet(CACHE/'matched-with-gaia-covariance.parquet')
roles=pd.read_parquet(CACHE/'stellar-spatial-holdouts.parquet')
d=d.merge(roles[['source_id','holdout_role','sky_pixel']],on='source_id',validate='one_to_one')
d=d[d.holdout_role=='training'].copy()
a=coordinates(*[d[k].to_numpy(float) for k in ['ra','dec','dist50','pmra','pmdec','VHELIO_AVG']])
d['R']=a[:,0];d['z']=a[:,1];d['phi']=(np.degrees(a[:,2])-27+180)%360-180
d['vR']=a[:,3];d['vphi']=a[:,4];d['vz']=a[:,5]
d['region']=region(d.R.to_numpy(),d.z.to_numpy())
d['Rbin']=np.floor(d.R/.5).astype(int)
d['azbin']=np.floor((d.phi+180)/45).astype(int)
d['FeHbin']=pd.cut(d.FE_H,[-3,-.5,0,1],include_lowest=True).astype(str)
d['alpha_bin']=(d.ALPHA_M>=.15).astype(int)
keys=['Rbin','azbin','FeHbin','alpha_bin']
out=[]
for pair in [(1,2),(3,4)]:
    sub=d[d.region.isin(pair)]
    bins=[]
    for key,g in sub.groupby(keys,observed=True):
        parts=[g[g.region==k] for k in pair]
        if min(len(p) for p in parts)==0:continue
        bins.append((key,parts))
    for minimum in [3,5,10]:
        common=[(key,parts) for key,parts in bins if min(len(p) for p in parts)>=minimum]
        if not common:
            out.append({'regions':pair,'minimum_per_region_cell':minimum,'common_cells':0});continue
        weights=np.array([min(len(p) for p in parts) for _,parts in common],float)
        weights/=weights.sum()
        sides=[]
        for side in [0,1]:
            mean=np.array([[parts[side][v].mean() for v in ['vR','vphi','vz']] for _,parts in common])
            second=np.array([[(parts[side][v]**2).mean() for v in ['vR','vphi','vz']] for _,parts in common])
            avg=np.sum(weights[:,None]*mean,axis=0)
            std=np.sqrt(np.maximum(0,np.sum(weights[:,None]*second,axis=0)-avg**2))
            ids=pd.concat([parts[side] for _,parts in common])
            sides.append({'region':pair[side],'stars_in_common_cells':len(ids),'sky_pixels':int(ids.sky_pixel.nunique()),
                          'weighted_means_kms':avg.tolist(),'weighted_stds_kms':std.tolist()})
        out.append({'regions':pair,'minimum_per_region_cell':minimum,'common_cells':len(common),
            'velocity_order':['vR','vphi','vz'],'sides':sides,
            'off_plane_minus_plane_mean_kms':(np.array(sides[1]['weighted_means_kms'])-sides[0]['weighted_means_kms']).tolist(),
            'status':'Descriptive training comparison standardized over coarse shared cells; no intrinsic-dispersion inference, selection correction, or confidence claim'})
save('training-overlap-comparison.json',out)
for x in out:print(x)
