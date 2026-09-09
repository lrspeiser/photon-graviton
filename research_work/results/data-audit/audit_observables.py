from pathlib import Path
import os
import json,csv,hashlib,gzip,collections,zipfile
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',Path(__file__).resolve().parents[3]/'research_work/generated'))/'data-audit'
OUT.mkdir(parents=True,exist_ok=True)
paths={
'sparc_catalog':'companion_wave_test/data/SPARC_Lelli2016c.mrt',
'sparc_curves':'companion_wave_test/data/Rotmod_LTG.zip',
'sparc_split':'companion_wave_test/data/sparc_frozen.json',
'dustpedia_themis':'companion_causal_test/data/themis.dat',
'dustpedia_dl14':'companion_causal_test/data/dl14.dat',
'dustpedia_schema':'companion_causal_test/data/ReadMe',
'matched_energy_derived':'companion_causal_test/energy.csv',
'cluster_derived':'companion_wave_test/cluster_comparison.csv',
'cluster_source_text':'companion_wave_test/data/xcop.txt',
'cluster_gas_source_text':'companion_wave_test/data/xcop_gas.txt',
'supernova_standardized':'shared_interaction_test/data/Pantheon+SH0ES.dat',
'supernova_covariance':'shared_interaction_test/data/Pantheon+SH0ES_STAT+SYS.cov.gz',
'supernova_schema':'shared_interaction_test/data/README',
'transient_widths':'time_revision/data/DES_event_averages.csv',
'spectral_aging':'time_first_principles/data/spectral_aging.json',
'cmb_spectrum':'time_first_principles/data/firas.txt',
'cmb_TT':'temporal_candidate_audit/data/planck_TT.txt',
'cmb_EE':'temporal_candidate_audit/data/planck_EE.txt',
'cf4_catalog':'temporal_candidate_audit/data/cf4_table2.dat',
'cf4_schema':'temporal_candidate_audit/data/cf4_ReadMe.txt',
'clock_paper_text':'nonrotating_action_test/data/atom_cavity.txt',
}
inventory=[]
for key,name in paths.items():
    p=ROOT/name
    row=dict(id=key,path=name,exists=p.is_file(),role='candidate input; inspect transformation notes before use')
    if p.is_file():row.update(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
    inventory.append(row)
catalog={}
for line in (ROOT/paths['sparc_catalog']).read_text(encoding='utf-8').splitlines():
    fields=line.split()
    if len(fields)!=19:continue
    try:
        method=int(fields[4]);dist=float(fields[2]);inc=float(fields[5]);quality=int(fields[17])
    except ValueError:continue
    name=fields[0]
    if method not in range(1,6):continue
    catalog[name]=dict(galaxy=name,distance_Mpc=dist,distance_method=method,inclination_deg=inc,quality=quality)
assert len(catalog)==175
split=json.loads((ROOT/paths['sparc_split']).read_text(encoding='utf-8'))['split']
selected=[n for names in split.values() for n in names]
assert len(selected)==len(set(selected))==149
matched=[r['galaxy'] for r in csv.DictReader((ROOT/paths['matched_energy_derived']).open(encoding='utf-8-sig'))]
assert len(matched)==len(set(matched))==26 and set(matched)<=set(selected)
methods={1:'Hubble flow with H0=73 and Virgo infall correction',2:'Tip of red giant branch magnitude',3:'Cepheid magnitude-period relation',4:'Ursa Major group distance',5:'Supernova light curve'}
counts={}
for label,names in [('catalog',list(catalog)),('analysis',selected),('matched',matched)]:
    count=collections.Counter(catalog[n]['distance_method'] for n in names)
    counts[label]={'total':len(names),'methods':{str(k):dict(count=count[k],description=v) for k,v in methods.items()}}
with zipfile.ZipFile(ROOT/paths['sparc_curves']) as archive:
    curves=[n for n in archive.namelist() if n.endswith('_rotmod.dat')]
    curve_header=archive.read(curves[0]).decode('utf-8').splitlines()[:4]
    assert len(curves)==175
manifest=json.loads((ROOT/'shared_interaction_test/data/manifest.json').read_text(encoding='utf-8'))
expected={m['name']:m['sha256'] for m in manifest}
dat=ROOT/paths['supernova_standardized'];cov=ROOT/paths['supernova_covariance']
covbytes=gzip.decompress(cov.read_bytes())
covlines=covbytes.splitlines();dimension=int(covlines[0]);entries=len(covlines)-1
actual_rows=len(dat.read_text(encoding='utf-8').splitlines())-1
assert entries==dimension**2 and actual_rows==dimension
assert hashlib.sha256(dat.read_bytes()).hexdigest()==expected[dat.name]
assert hashlib.sha256(covbytes).hexdigest()==expected['Pantheon+SH0ES_STAT+SYS.cov']
result=dict(scope='Selected locally recovered input products; not an exhaustive raw-observation archive or full original-manifest re-audit',inventory=inventory,distance_method_counts=counts,matched_distance_details=[catalog[n] for n in matched],sparc_curve_members=len(curves),curve_header=curve_header,pantheon=dict(rows=actual_rows,covariance_dimension=dimension,covariance_entries=entries,data_and_decompressed_covariance_match_archived_hashes=True),all_listed_paths_present=all(r['exists'] for r in inventory))
(OUT/'input-audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps({k:result[k] for k in ['distance_method_counts','pantheon','all_listed_paths_present']},indent=2))
