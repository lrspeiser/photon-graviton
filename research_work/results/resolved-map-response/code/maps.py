"""Measured map extraction and positive radial/angular source reconstruction."""
from __future__ import annotations
import json, argparse, hashlib
from pathlib import Path
import numpy as np
from scipy.linalg import solve
from fits_images import read_images,channel

IDS=['7443-12703','8082-12704','8081-1901']

def basis(r,phi,limit,mmax=2):
    c=np.linspace(0,limit,6);bw=limit/4
    rb=np.exp(-.5*((np.asarray(r)[...,None]-c)/bw)**2);rb/=rb.sum(-1,keepdims=True)
    terms=[np.ones_like(phi)]
    for m in range(1,mmax+1):terms += [np.cos(m*phi),np.sin(m*phi)]
    return (rb[..., :,None]*np.stack(terms,-1)[...,None,:]).reshape(np.shape(r)+(-1,))

def log_fit(r,phi,flux,snr,limit):
    X=basis(r,phi,limit);w=np.minimum(snr,20.)**2;w/=w.mean()
    # Fixed mild ridge on radial/azimuthal coefficients; no outcome selection.
    penalty=np.tile(np.array([.005,.10,.10,.40,.40]),6)
    normal=X.T@(w[:,None]*X)+np.diag(penalty)
    b=np.linalg.solve(normal,X.T@(w*np.log(flux)))
    return b, float(np.sqrt(np.mean((X@b-np.log(flux))**2)))

def shape(coeff,r,phi,limit):
    return np.exp(basis(np.asarray(r),np.asarray(phi),limit)@coeff)

def velocity_fit(h,common,limit):
    rc=h['SPX_ELLCOO'][1][1];theta=h['SPX_ELLCOO'][1][3]*np.pi/180
    sky=h['SPX_SKYCOO'][1];bins=h['BINID'][1][1]
    ha=channel(h['EMLINE_GVEL'][0],'Ha-6564')
    sv=h['STELLAR_VEL'][1];gv=h['EMLINE_GVEL'][1][ha]
    si=h['STELLAR_VEL_IVAR'][1];gi=h['EMLINE_GVEL_IVAR'][1][ha]
    valid=common&(h['STELLAR_VEL_MASK'][1]==0)&(h['EMLINE_GVEL_MASK'][1][ha]==0)&(si>0)&(gi>0)&(bins>=0)&(rc<=limit)
    # Flexible first-harmonic velocity field. Not a gravitational model.
    centers=np.linspace(0,limit,5);width=limit/3
    rb=np.exp(-.5*((rc[...,None]-centers)/width)**2);rb/=rb.sum(-1,keepdims=True)
    design=np.concatenate([np.ones(rc.shape+(1,)),rb*np.cos(theta)[...,None],rb*np.sin(theta)[...,None]],axis=-1)
    x=[];ys=[];yg=[];ws=[];wg=[];xy=[];bin_ids=[]
    light=np.maximum(h['SPX_MFLUX'][1],0)
    for b in np.unique(bins[valid]):
        mask=valid&(bins==b); ww=light[mask];ww/=ww.sum()
        x.append(ww@design[mask]);ys.append(np.median(sv[mask]));yg.append(ww@gv[mask]);ws.append(np.median(si[mask]));wg.append(np.median(gi[mask]));xy.append([1.,ww@sky[0][mask],ww@sky[1][mask]]);bin_ids.append(int(b))
    X=np.asarray(x);xy=np.asarray(xy);ys=np.asarray(ys);yg=np.asarray(yg)
    out={'valid_spaxels':int(valid.sum()),'unique_stellar_bins':len(ys),'scope':'Binned projected velocity harmonics; gas bin means use median inverse variance as a working weight, not independent-pixel covariance'}
    co=[]
    for name,y,w in [('star',ys,np.array(ws)),('gas',yg,np.array(wg))]:
        w=w/np.median(w);reg=np.eye(X.shape[1])*.01;reg[0,0]=0
        coef=np.linalg.solve(X.T@(w[:,None]*X)+reg,X.T@(w*y));plane=np.linalg.lstsq(xy*np.sqrt(w[:,None]),y*np.sqrt(w),rcond=None)[0]
        out[name]={'harmonic_coefficients':coef.tolist(),'plane_coefficients':plane.tolist(),'plane_RMS_kms':float(np.sqrt(np.mean((xy@plane-y)**2))),'harmonic_RMS_kms':float(np.sqrt(np.mean((X@coef-y)**2)))};co.append(coef)
    sample=np.array([.5,.8,1.1,1.4])*limit/1.5
    rbas=np.exp(-.5*((sample[:,None]-centers)/width)**2);rbas/=rbas.sum(1,keepdims=True)
    amp=[];phase=[]
    for coef in co:
        a=rbas@coef[1:6];b=rbas@coef[6:];amp.append(np.hypot(a,b));phase.append(np.arctan2(b,a))
    delta=np.angle(np.exp(1j*(phase[1]-phase[0])))
    out['radial_samples']=[dict(R_over_Re=float(r),star_first_harmonic_kms=float(a),gas_first_harmonic_kms=float(b),projected_phase_difference_deg=float(d*180/np.pi)) for r,a,b,d in zip(sample,*amp,delta)]
    return out,dict(bin_design=X,bin_sky_design=xy,star_bin_velocity=ys,gas_bin_velocity=yg,bin_ids=np.array(bin_ids),star_fit=X@co[0],gas_fit=X@co[1])

def prepare(inputs,out):
    out.mkdir(parents=True,exist_ok=False);records=[]
    for ident in IDS:
        f=next(inputs.glob(f'manga-{ident}-*.fits.gz'));h=read_images(f);hdr=h['PRIMARY'][0]
        old=np.load(inputs/f'{ident}-maps.npz');ha=channel(h['EMLINE_GFLUX'][0],'Ha-6564')
        for current,key in [(h['STELLAR_VEL'][1],'star_velocity'),(h['EMLINE_GFLUX'][1][ha],'Halpha'),(h['SPX_MFLUX'][1],'light')]:
            assert np.array_equal(current,old[key]),key
        L=h['SPX_MFLUX'][1].astype(float);H=h['EMLINE_GFLUX'][1][ha].astype(float)
        hs=H*np.sqrt(np.maximum(h['EMLINE_GFLUX_IVAR'][1][ha],0));ls=h['SPX_SNR'][1]
        common=(L>0)&(H>0)&(hs>=5)&(ls>=5)&(h['EMLINE_GFLUX_MASK'][1][ha]==0)
        r=h['SPX_ELLCOO'][1][1];phi=np.deg2rad(h['SPX_ELLCOO'][1][3]);limit=float(min(2.,np.percentile(r[common],90)))
        mask=common&(r<=limit)
        cb,cl=log_fit(r[mask],phi[mask],L[mask],ls[mask],limit);gb,gl=log_fit(r[mask],phi[mask],H[mask],hs[mask],limit)
        velocity,arrays=velocity_fit(h,common,limit)
        # All spectral ratios stay observed line ratios, not temperatures/densities.
        ratios={}
        for a,b in [('Ha-6564','Hb-4862'),('NII-6585','Ha-6564'),('SII-6718','SII-6732')]:
            ia,ib=[channel(h['EMLINE_GFLUX'][0],s) for s in (a,b)];fa,fb=h['EMLINE_GFLUX'][1][[ia,ib]]
            ma=(h['EMLINE_GFLUX_MASK'][1][ia]==0)&(h['EMLINE_GFLUX_MASK'][1][ib]==0)&mask
            ma &= (fa*np.sqrt(np.maximum(h['EMLINE_GFLUX_IVAR'][1][ia],0))>=5)&(fb*np.sqrt(np.maximum(h['EMLINE_GFLUX_IVAR'][1][ib],0))>=5)&(fb>0)
            ratios[a+'/'+b]=dict(pixels=int(ma.sum()),median=float(np.median(fa[ma]/fb[ma])) if ma.any() else None)
        # Counterfactual relative rotations of positive fitted projected shapes.
        rr=(np.arange(60)+.5)*limit/60;pp=(np.arange(144)+.5)*2*np.pi/144
        R,F=np.meshgrid(rr,pp,indexing='ij');lp=shape(cb,R,F,limit);hp=shape(gb,R,F,limit)
        lres=lp/lp.mean(1,keepdims=True)-1;hres=hp/hp.mean(1,keepdims=True)-1;w=R/R.sum()
        overlaps=[float(np.sum(w*lres*np.roll(hres,k,axis=1))) for k in range(144)]
        curves=dict(phi_shift_degrees=(np.arange(144)*2.5).tolist(),radially_detrended_overlap=overlaps,observed=overlaps[0],rotation_percentile=float(np.mean(np.asarray(overlaps)<=overlaps[0])),scope='Descriptive phase-control percentile, not a probability under an observational null')
        rec=dict(id=ident,input_sha256=hashlib.sha256(f.read_bytes()).hexdigest(),map_shape=list(L.shape),valid_photometry_pixels=int(mask.sum()),aperture_Re=limit,Re_arcsec=hdr['REFF'],assumed_axis_ratio=1-hdr['ECOOELL'],photometric_PA_deg=hdr['ECOOPA'],beam_FWHM_arcsec=max(hdr['GFWHM'],hdr['RFWHM']),beam_FWHM_over_Re=max(hdr['GFWHM'],hdr['RFWHM'])/hdr['REFF'],source_log_fit_RMS=cl,gas_proxy_log_fit_RMS=gl,source_coeff=cb.tolist(),gas_coeff=gb.tolist(),velocity=velocity,observed_line_ratios=ratios,overlap=curves)
        records.append(rec)
        np.savez_compressed(out/(ident+'.npz'),light=L,Halpha=H,mask=mask,r_Re=r,phi=phi,sky=h['SPX_SKYCOO'][1],source_coeff=cb,gas_coeff=gb,**arrays)
        print(ident,'bins',velocity['unique_stellar_bins'],'overlap',curves['observed'],'rank',curves['rotation_percentile'],'kin',velocity['radial_samples'],flush=True)
    (out/'maps.json').write_text(json.dumps(records,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':
    ap=argparse.ArgumentParser(__doc__);ap.add_argument('--inputs',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();prepare(a.inputs,a.output)
