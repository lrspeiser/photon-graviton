"""Acquire source manifest and reconstruct published blank-field controls."""
from pathlib import Path
import hashlib,json,tarfile
import requests,fitz
from scipy.stats import chi2

OUT=Path(__file__).resolve().parent
CACHE=OUT.parents[2]/'research_work/generated/coma-kubo2007'
archive=CACHE/'arxiv-source'
url='https://arxiv.org/src/0709.0506'
if not archive.exists():
    r=requests.get(url,timeout=30);r.raise_for_status();archive.write_bytes(r.content)
manifest=[]
with tarfile.open(archive) as tar:
    for member in tar.getmembers():
        assert member.isfile()
        data=tar.extractfile(member).read()
        manifest.append(dict(name=member.name,bytes=len(data),sha256=hashlib.sha256(data).hexdigest()))
assert {m['name'] for m in manifest}=={'ms.tex','f1.eps','f2.eps','f3.eps'}
prior=json.loads((OUT/'kubo-figure-data.json').read_text())
pdf=CACHE/'paper.pdf'
assert hashlib.sha256(pdf.read_bytes()).hexdigest()==prior['pdf_sha256']
page=fitz.open(pdf)[15]
assert 'Fig. 3.' in page.get_text()
drawings=page.get_drawings()
assert len(drawings)==8
lines=drawings[7]['items'];assert len(lines)==67
baselines=[line for line in lines if line[0]=='l' and abs(line[1].y-line[2].y)<.001 and line[2].x-line[1].x>400]
assert len(baselines)==1
x0,y0=baselines[0][1];right=baselines[0][2].x
bottom=max(p.y for p in drawings[0]['items'][0][1])
xs=(right-x0)/10.5;ys=(bottom-y0)/.005
assert 39<xs<40 and 24000<ys<25000
rows=[]
for i in range(6):
    up,lo=lines[2*i],lines[12+2*i]
    cross_up,cross_lo=lines[43+2*i],lines[55+2*i]
    x,y=up[1]
    assert abs(cross_up[1].x-x)<.001 and abs(cross_lo[1].x-x)<.001
    assert cross_up[2].y<cross_up[1].y<cross_lo[2].y
    marker=drawings[i+1]['rect']
    assert abs((marker.x0+marker.x1)/2-x)<.21 and abs((marker.y0+marker.y1)/2-y)<.21
    rows.append(dict(bin=i+1,published_radius_h_inverse_Mpc=(x-x0)/xs,
        shear_t=(y0-y)/ys,sigma_t=(lo[2].y-up[2].y)/(2*ys),
        shear_cross=(y0-cross_up[1].y)/ys,sigma_cross=(cross_lo[2].y-cross_up[2].y)/(2*ys)))
ct=sum((r['shear_t']/r['sigma_t'])**2 for r in rows)
cx=sum((r['shear_cross']/r['sigma_cross'])**2 for r in rows)
assert abs(ct/1.84-1)<.02 and abs(cx/12.58-1)<.02
# Conditional offset subtraction, retaining uncertainty in the control estimate.
adjusted=[]
for source,blank in zip(prior['rows'],rows):
    # Same bin ID, not identical representative radius. No interpolation is implied.
    var=source['plotted_sigma_t']**2+blank['sigma_t']**2
    adjusted.append(dict(bin=source['bin'],shear_t=source['shear_t']-blank['shear_t'],sigma=var**.5))
adjusted_chi=sum((r['shear_t']/r['sigma'])**2 for r in adjusted)
result=dict(source_archive_url=url,archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),archive_members=manifest,
    source_audit='Archive contains manuscript and three EPS figures only; no separate shear table, pipeline or covariance file. No exact intermediate angular bin boundaries recovered.',
    blank_figure=3,page_one_based=16,pdf_sha256=prior['pdf_sha256'],blank_rows=rows,
    checks=dict(blank_t_null_chi2=ct,paper_blank_t_null_chi2=1.84,blank_cross_null_chi2=cx,paper_blank_cross_null_chi2=12.58),
    conditional_offset_subtraction=dict(assumptions='Common additive offset by bin estimated from published blank mean; independent Coma/control errors; diagonal bins. Sensitivity analysis only, not validated covariance or required correction.',
        rows=adjusted,null_chi2=adjusted_chi,degrees_of_freedom=6,conditional_null_p=float(chi2.sf(adjusted_chi,6))),
    status='Exploratory reconstructed control means. Six individual blank-field profiles remain unavailable; do not infer full covariance from their published mean.')
(OUT/'kubo-controls-results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(dict(checks=result['checks'],conditional=result['conditional_offset_subtraction']),indent=2))
