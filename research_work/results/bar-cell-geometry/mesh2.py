"""Second nested geodesic subdivision; new directions need actual trajectories."""
from pathlib import Path
import json
import numpy as np
from scipy.spatial import ConvexHull
HERE=Path(__file__).resolve().parent
SYM=np.array([[1,1,1],[-1,-1,1],[1,1,-1],[-1,-1,-1]])

def canonical(direction):
    a=np.array(direction,dtype=float);a[abs(a)<1e-14]=0.;a/=np.linalg.norm(a)
    for s in SYM:
        b=a*s
        if b[2]>=0 and (b[1]>0 or (b[1]==0 and b[0]>=0)):
            return b,s
    raise ValueError('No symmetry representative')

def key(direction):return tuple(np.round(direction,12))

if __name__=='__main__':
    source=json.loads((HERE/'mesh.json').read_text(encoding='utf8'))
    vertices=source['vertices'].copy()
    parent=np.array(source['faces'])
    edge_indices={};children=[]
    for a,b,c in parent:
        mids=[]
        for i,j in ((a,b),(b,c),(a,c)):
            edge=tuple(sorted((int(i),int(j))))
            if edge not in edge_indices:
                v=np.array(vertices[i])+vertices[j];v/=np.linalg.norm(v)
                edge_indices[edge]=len(vertices);vertices.append(v.tolist())
            mids.append(edge_indices[edge])
        ab,bc,ac=mids
        children.extend([[int(a),ab,ac],[ab,int(b),bc],[ac,bc,int(c)],[ab,bc,ac]])
    vertices=np.array(vertices);faces=np.array(children)
    a,b,c=np.moveaxis(vertices[faces],1,0)
    area=2*np.arctan2(abs(np.einsum('ij,ij->i',a,np.cross(b,c))),1+np.sum(a*b+b*c+c*a,axis=1))
    assert abs(sum(area)-4*np.pi)<1e-10
    unique={};unique_values=[];lookup=[];signs=[]
    for vertex in vertices:
        v,s=canonical(vertex);k=key(v)
        if k not in unique:
            unique[k]=len(unique);unique_values.append(v.tolist())
        lookup.append(unique[k]);signs.append(s.tolist())
    directions=unique_values
    # Keys are rounded only for matching; actual direction coordinates are retained.
    directions=np.array(directions);directions/=np.linalg.norm(directions,axis=1)[:,None]
    error=float(np.max(np.linalg.norm(directions[lookup]*signs-vertices,axis=1)))
    assert error<1e-12
    out=dict(vertices=vertices.tolist(),faces=faces.tolist(),parent_faces=parent.tolist(),
             directions=directions.tolist(),lookup=lookup,signs=signs,reconstruction_error=error,
             solid_angle=float(sum(area)),parent_vertices=len(source['vertices']))
    (HERE/'mesh2.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
    print('Vertices',len(vertices),'faces',len(faces),'unique trajectories per source',len(directions),'error',error)
