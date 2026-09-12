"""Conforming local angular refinement; retain every unselected source region."""
from pathlib import Path
import json,hashlib
import numpy as np
from mesh import canonical,key
H=Path(__file__).resolve().parent
source=json.loads((H/'mesh2.json').read_text());ranking=json.loads((H/'angular-attribution-summary.json').read_text())
selected_parent={r['parent_face'] for r in ranking['ranked_regions'][:ranking['counts_for_norm_fraction']['0.9']]}
base=np.array(source['faces']);vertices=[list(x) for x in source['vertices']];edges={}
for i,face in enumerate(base):
 if i//4 not in selected_parent:continue
 for a,b in zip(face,np.roll(face,-1)):
  e=tuple(sorted((int(a),int(b))))
  if e not in edges:
   v=np.array(vertices[a])+vertices[b];v/=np.linalg.norm(v);edges[e]=len(vertices);vertices.append(v.tolist())
children=[];parents=[]
for i,face in enumerate(base):
 v=list(map(int,face));mid=[edges.get(tuple(sorted((v[k],v[(k+1)%3])))) for k in range(3)]
 count=sum(m is not None for m in mid)
 if count==0:tri=[v]
 elif count==3:
  a,b,c=v;ab,bc,ca=mid;tri=[[a,ab,ca],[ab,b,bc],[ca,bc,c],[ab,bc,ca]]
 elif count==1:
  k=next(k for k in range(3) if mid[k] is not None);a,b,c=[v[(k+j)%3] for j in range(3)];m=mid[k]
  tri=[[a,m,c],[m,b,c]]
 else:
  k=next(k for k in range(3) if mid[k] is not None and mid[(k+1)%3] is not None)
  a,b,c=[v[(k+j)%3] for j in range(3)];ab=mid[k];bc=mid[(k+1)%3]
  tri=[[b,bc,ab],[a,ab,c],[ab,bc,c]]
 children.extend(tri);parents.extend([i]*len(tri))
V=np.array(vertices);F=np.array(children);P=np.array(parents)
a,b,c=np.moveaxis(V[F],1,0)
omega=2*np.arctan2(abs(np.einsum('ij,ij->i',a,np.cross(b,c))),1+np.sum(a*b+b*c+c*a,axis=1))
a,b,c=np.moveaxis(V[base],1,0)
oldomega=2*np.arctan2(abs(np.einsum('ij,ij->i',a,np.cross(b,c))),1+np.sum(a*b+b*c+c*a,axis=1))
areaerror=float(np.max(abs(np.bincount(P,weights=omega,minlength=len(base))-oldomega)))
assert areaerror<1e-12 and abs(omega.sum()-4*np.pi)<1e-10
E=np.sort(np.concatenate([F[:,[0,1]],F[:,[1,2]],F[:,[2,0]]]),axis=1);unique_edges,counts=np.unique(E,axis=0,return_counts=True)
assert np.all(counts==2) and len(V)-len(unique_edges)+len(F)==2
unique={};values=[];lookup=[];signs=[]
for vertex in V:
 v,sign=canonical(vertex);k=key(v)
 if k not in unique:unique[k]=len(values);values.append(v.tolist())
 lookup.append(unique[k]);signs.append(sign.tolist())
error=float(np.max(np.linalg.norm(np.array(values)[lookup]*signs-V,axis=1)));assert error<1e-12
out=dict(vertices=vertices,faces=children,directions=values,lookup=lookup,signs=signs,parent_second=parents,
 selected_first_parents=sorted(selected_parent),parent_hash=hashlib.sha256((H/'mesh2.json').read_bytes()).hexdigest(),
 checks=dict(parent_solid_angle_error=areaerror,solid_angle=float(omega.sum()),euler=2,edges_have_two_faces=True,reconstruction_error=error))
(H/'mesh-adaptive.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
print('vertices',len(V),'faces',len(F),'representatives',len(values),'new',len(values)-len(source['directions']),out['checks'])
