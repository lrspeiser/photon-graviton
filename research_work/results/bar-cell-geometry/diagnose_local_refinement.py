"""Locate resolved interpolation errors; distinguish unsampled cells from passing cells."""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.interpolate import CubicHermiteSpline
H=Path(__file__).resolve().parent;ROOT=H.parents[2]
M=json.loads((H/'mesh.json').read_text()); V=np.array(M['vertices']); F=np.array(M['faces'])
L=np.array(M['lookup']); S=np.array(M['signs']); SYM=np.array([[1,1,1],[-1,-1,1],[1,1,-1],[-1,-1,-1]])
A=np.linalg.inv(V[F].transpose(0,2,1)); times=np.array([0,.005,.01,.025,.05,.1,.25])
a,b,c=np.moveaxis(V[F],1,0)
omega=2*np.arctan2(abs(np.einsum('ij,ij->i',a,np.cross(b,c))),1+np.sum(a*b+b*c+c*a,axis=1))
old=json.loads((H.parent/'bar-volume-gravity/prepared8.json').read_text())
def spline(r):
 p=ROOT/r['cache'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['cache_sha256']
 d=np.load(p);return CubicHermiteSpline(d['ages'],d['position'],d['velocity'])(times)
out=[]
for R in (1,3):
 new=json.loads((H/f'prepared-R{R}.json').read_text())['records']
 C=np.stack([spline(r) for r in new],axis=1)[:,L]*S[None,:,:]
 cap=np.array([r['capture'] for r in new])[L]
 mass=omega*np.mean(cap[F],axis=1);mass/=mass.sum()
 selected=[];errors=[];weights=[]
 for r in [r for r in old['records'] if r['R']==R]:
  actual=spline(r)
  for sign in SYM:
   direction=np.array(r['direction'])*sign
   alpha=np.einsum('ijk,k->ij',A,direction)
   ids=np.flatnonzero(np.all(alpha>=-1e-10,axis=1)&(alpha.sum(axis=1)>0));assert len(ids)
   i=int(ids[0]);w=alpha[i]/alpha[i].sum();assert abs(w.sum()-1)<1e-12
   predicted=np.sum(C[:,F[i]]*w[None,:,None],axis=1)
   selected.append(i);errors.append(np.linalg.norm(predicted-actual*sign,axis=1));weights.append(r['raw_weight']/4)
 errors=np.array(errors).T;weights=np.array(weights);weights/=weights.sum()
 occupied=np.unique(selected);records=[]
 for j,T in enumerate(times):
  contribution=weights*errors[j]**2
  cell=np.bincount(selected,weights=contribution,minlength=len(F))
  order=np.argsort(-cell);total=cell.sum()
  counts={str(frac):int(np.searchsorted(np.cumsum(cell[order]),frac*total)+1) for frac in (.5,.9,.99)}
  edges=np.stack([C[j,F[:,1]]-C[j,F[:,0]],C[j,F[:,2]]-C[j,F[:,0]],C[j,F[:,2]]-C[j,F[:,1]]],axis=1)
  extent=np.max(np.linalg.norm(edges,axis=2),axis=1)/R
  records.append(dict(T=float(T),rms_over_R=float(np.sqrt(total)/R),cells_for_error_fraction=counts,
   sampled_weight_above_002R=float(weights[errors[j]>.02*R].sum()),
   source_mass_with_cell_diameter_above_R=float(mass[extent>1].sum()),
   ranked_cells=[dict(face=int(i),squared_error_fraction=float(cell[i]/total),source_mass=float(mass[i])) for i in order if cell[i]>0]))
 out.append(dict(R=R,faces=len(F),reference_directions=len(selected),sampled_faces=len(occupied),
  sampled_face_source_mass=float(mass[occupied].sum()),records=records))
(H/'local-refinement-diagnostic.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
for r in out:
 print({k:v for k,v in r.items() if k!='records'})
 print([{k:v for k,v in t.items() if k!='ranked_cells'} for t in r['records']])
