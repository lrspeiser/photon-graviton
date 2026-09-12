from pathlib import Path
import json,csv
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad,cumulative_trapezoid

HERE=Path(__file__).resolve().parent
rows=[]
for R in [10.,30.]:
    for count,order in [(321,48),(641,96)]:
        r=np.linspace(0,R,count);mu,wm=leggauss(order);z,w=leggauss(order);u=(z+1)/2;wu=w/2
        rr=r[:,None,None];mm=mu[None,:,None]
        sb=rr*mm+np.sqrt(R*R-rr*rr+rr*rr*mm*mm);s=sb*u[None,None,:]
        radii=np.sqrt(np.maximum(0,rr*rr+s*s-2*rr*mm*s))
        impact=R*u;half=np.sqrt(R*R-impact*impact)
        chord_r=np.sqrt(impact[:,None]**2+(half[:,None]*u[None,:])**2)
        for p in [1.,2.,3.]:
            shape=lambda r:(2/(1+np.sqrt(1+r*r)))**p
            half_depth=quad(shape,0,R,epsabs=1e-12)[0]
            base=sb[:,:,0]*np.sum(wu*shape(radii),axis=2)
            chord_base=2*half*np.sum(wu*shape(chord_r),axis=1)
            for depth in [.1,1.,10.]:
                A=depth/(2*half_depth);kappa=A*shape(r)
                intensity=np.sum(np.exp(-A*base)*wm[None,:]/2,axis=1)
                center_error=abs(intensity[0]/np.exp(-depth/2)-1)
                q=4*np.pi*kappa*intensity
                enc=4*np.pi*cumulative_trapezoid(q*r*r,r,initial=0)
                absorbed=8*np.pi*np.pi*R*np.sum(wu*impact*(-np.expm1(-A*chord_base)))
                power_error=abs(enc[-1]/absorbed-1)
                if count==641:assert center_error<1e-6 and power_error<.002
                points=np.linspace(3,9,121);masses=np.interp(points,r,enc)
                v=np.sqrt(masses/points)
                slope=float(np.polyfit(np.log(points),np.log(v),1)[0])
                density_slope=float(np.polyfit(np.log(points),np.log(np.interp(points,r,q)),1)[0])
                b=5.
                def piecewise(fn,edges,n):
                    nodes,weights=leggauss(n);left=edges[:-1,None];width=np.diff(edges)[:,None]
                    return float(np.sum(fn(left+width*(nodes+1)/2)*weights*width/2))
                inner_edges=np.unique(np.r_[0,r[(r>0)&(r<b)],b]);outer_edges=np.sqrt(np.unique(np.r_[b,r[r>b]])**2-b*b)
                inner=piecewise(lambda t:4*np.pi*t*t*np.interp(t,r,q),inner_edges,4)
                def outer_fn(t):
                    s=np.sqrt(b*b+t*t);return 4*np.pi*np.interp(s,r,q)*t*(s-t)
                outer=piecewise(outer_fn,outer_edges,8);check=piecewise(outer_fn,outer_edges,16)
                quaderr=abs(check-outer)/(inner+outer);assert quaderr<1e-9
                rows.append(dict(outer_radius=R,radial_nodes=count,quadrature_order=order,p=p,diameter_optical_depth=depth,kappa_center=A,
                    captured_fraction=float(absorbed/(4*np.pi*np.pi*R*R)),power_agreement=float(power_error),central_intensity_error=float(center_error),
                    v9_over_v3=float(v[-1]/v[0]),log_speed_slope=slope,log_density_slope=density_slope,
                    v27_over_v9=float(np.sqrt((np.interp(27,r,enc)/27)/(np.interp(9,r,enc)/9))) if R==30 else None,
                    enclosed_fraction_3=float(np.interp(3,r,enc)/enc[-1]),enclosed_fraction_5=float(np.interp(5,r,enc)/enc[-1]),enclosed_fraction_9=float(np.interp(9,r,enc)/enc[-1]),
                    projected_over_spherical_mass_at_5=float((inner+outer)/np.interp(5,r,enc)),projection_quadrature_error=float(quaderr)))
ref=[];boundary=[]
for R in [10.,30.]:
    for p in [1.,2.,3.]:
        for depth in [.1,1.,10.]:
            lo=next(a for a in rows if a['outer_radius']==R and a['p']==p and a['diameter_optical_depth']==depth and a['radial_nodes']==321)
            hi=next(a for a in rows if a['outer_radius']==R and a['p']==p and a['diameter_optical_depth']==depth and a['radial_nodes']==641)
            dv=abs(lo['v9_over_v3']-hi['v9_over_v3']);ds=abs(lo['log_speed_slope']-hi['log_speed_slope'])
            assert dv<.003 and ds<.003
            ref.append(dict(R=R,p=p,depth=depth,v_ratio_difference=dv,slope_difference=ds))
for p in [1.,2.,3.]:
    for depth in [.1,1.,10.]:
        a=next(v for v in rows if v['outer_radius']==10 and v['p']==p and v['diameter_optical_depth']==depth and v['radial_nodes']==641)
        b=next(v for v in rows if v['outer_radius']==30 and v['p']==p and v['diameter_optical_depth']==depth and v['radial_nodes']==641)
        boundary.append(dict(p=p,depth=depth,v_ratio_change=b['v9_over_v3']-a['v9_over_v3'],lensing_ratio_change=b['projected_over_spherical_mass_at_5']-a['projected_over_spherical_mass_at_5']))
with (HERE/'cases.csv').open('w',encoding='utf8',newline='') as f:
    wr=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');wr.writeheader();wr.writerows(rows)
(HERE/'results.json').write_text(json.dumps(dict(cases=rows,refinement=ref,boundary_sensitivity=boundary,scope='Declared common opacity shapes; no astronomical fit or radiation-supply inference.'),indent=2)+'\n',encoding='utf8',newline='\n')
print(json.dumps([r for r in rows if r['radial_nodes']==641 and r['outer_radius']==10],indent=2))
