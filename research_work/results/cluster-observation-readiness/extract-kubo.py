"""Recover published figure coordinates, not a raw shear catalog or halo fit."""
from pathlib import Path
import json,hashlib
import requests,fitz

OUT=Path(__file__).resolve().parent
CACHE=OUT.parents[2]/'research_work/generated/coma-kubo2007'
CACHE.mkdir(parents=True,exist_ok=True)
url='https://lss.fnal.gov/archive/2007/pub/fermilab-pub-07-453-a-cd.pdf'
pdf=CACHE/'paper.pdf'
if not pdf.exists():
    response=requests.get(url,timeout=30);response.raise_for_status();pdf.write_bytes(response.content)
assert hashlib.sha256(pdf.read_bytes()).hexdigest()=='87f5126c3e64a03f37315623b811c20b8c7a37a8de40b46a9fb48e87e06f0f6a'
doc=fitz.open(pdf)
page=doc[14]
assert 'Fig. 2.' in page.get_text()
drawings=page.get_drawings()
assert len(drawings)==9 and len(drawings[7]['items'])==67
lines=drawings[7]['items']
baseline=lines[66]
x0,y0=baseline[1];xright=baseline[2].x
# PDF coordinates are top-down. Quad ordering is inspected rather than assumed.
ybottom=max(p.y for p in drawings[0]['items'][0][1])
# Visually verified axes: x left=0, right=10.5; lower y=-0.005; long line y=0.
xscale=(xright-x0)/10.5
yscale=(ybottom-y0)/.005
assert 39<xscale<40 and 24000<yscale<25000
rows=[]
for i in range(6):
    upper=lines[2*i];lower=lines[12+2*i]
    cupper=lines[42+2*i];clower=lines[54+2*i]
    assert upper[0]==lower[0]==cupper[0]==clower[0]=='l'
    x,y=upper[1]
    assert abs(lower[1].y-y)<.001 and abs(cupper[1].x-x)<.001
    marker=drawings[i+1]['rect']
    assert abs((marker.x0+marker.x1)/2-x)<.21 and abs((marker.y0+marker.y1)/2-y)<.21
    error=(lower[2].y-upper[2].y)/(2*yscale)
    cross_error=(clower[2].y-cupper[2].y)/(2*yscale)
    assert abs((y-upper[2].y)-(lower[2].y-y))<.21
    rows.append(dict(bin=i+1,published_radius_h_inverse_Mpc=(x-x0)/xscale,
        shear_t=(y0-y)/yscale,plotted_sigma_t=error,
        shear_cross=(y0-cupper[1].y)/yscale,plotted_sigma_cross=cross_error,
        pdf_center_points=[x,y],pdf_error_endpoints_y=[upper[2].y,lower[2].y]))
chi_t=sum((r['shear_t']/r['plotted_sigma_t'])**2 for r in rows)
chi_cross=sum((r['shear_cross']/r['plotted_sigma_cross'])**2 for r in rows)
assert abs(chi_t/23.33-1)<.01 and abs(chi_cross/5.65-1)<.01
result=dict(source=url,pdf_sha256=hashlib.sha256(pdf.read_bytes()).hexdigest(),page_one_based=15,figure=2,
    status='Exploratory vector-figure reconstruction; not author machine-readable measurements, raw shapes, full covariance or untouched validation.',
    axis_calibration=dict(x0=x0,xright=xright,x_units_per_point=1/xscale,y_zero=y0,y_bottom=ybottom,shear_per_point=1/yscale),
    conservative_plot_coordinate_tolerance_points=.25,plot_coordinate_tolerance_shear=.25/yscale,
    checks=dict(null_chi2_t=chi_t,paper_null_chi2_t=23.33,null_chi2_cross=chi_cross,paper_null_chi2_cross=5.65),
    restrictions=['Radius retains paper cosmological units; not adopted as a fictional physical distance.',
    'Do not substitute a fitted NFW curve for data; curve excluded from extraction.',
    'No bin edges, source weights, full covariance or source geometry recovered.',
    'Cross component and negative measured shear retained; do not clip.'],rows=rows)
(OUT/'kubo-figure-data.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
page.get_pixmap(matrix=fitz.Matrix(1.5,1.5)).save(CACHE/'figure-page.png')
print(json.dumps(result,indent=2))
