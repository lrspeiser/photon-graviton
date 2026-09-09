from pathlib import Path
import os
import json,hashlib,shutil,subprocess,sys,time,os,platform,importlib.metadata,csv,math,gzip
root=Path(__file__).resolve().parents[3]
base=root/'research_work';out=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',base/'generated'))/'baseline';out.mkdir(parents=True,exist_ok=True)
run=out/'scratch'
if (out/'run-status.json').exists():raise RuntimeError('Existing live or completed run: inspect its handle before restarting')
run.mkdir(parents=True,exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest=json.loads((root/'SNAPSHOT_MANIFEST.json').read_text(encoding='utf-8'));audit=[]
for row in manifest['files']:
 p=root/row['path'];compressed=p.with_name(p.name+'.gz');actual=sha(p) if p.exists() else hashlib.sha256(gzip.decompress(compressed.read_bytes())).hexdigest() if compressed.exists() else None
 audit.append({**row,'actual_sha256':actual,'status':'match' if actual==row['sha256'] else 'missing' if actual is None else 'mismatch'})
(out/'manifest-audit.json').write_text(json.dumps(audit,indent=2))
inputs=json.loads((root/'companion_causal_test/hashes.json').read_text(encoding='utf-8'))
assert all(sha(root/p)==h for p,h in inputs.items())
for folder in ['companion_causal_test','companion_wave_test','companion_deposition_fit']:shutil.copytree(root/folder,run/folder,dirs_exist_ok=True)
original_hashes={p.relative_to(root).as_posix():sha(p) for f in ['companion_causal_test','companion_wave_test','companion_deposition_fit'] for p in (root/f).rglob('*') if p.is_file()}
packages={}
for name in ['numpy','scipy','pandas','matplotlib','python-docx']:
 try: packages[name]=importlib.metadata.version(name)
 except importlib.metadata.PackageNotFoundError: packages[name]='not installed; not used by these numerical scripts'
(out/'environment-lock.txt').write_text('Python '+sys.version+'\nPlatform '+platform.platform()+'\n'+''.join(k+'=='+v+'\n' for k,v in packages.items())+'OPENBLAS_NUM_THREADS=1\nOMP_NUM_THREADS=1\nPYTHONUTF8=1\n')
# Declared before execution: numerical regression, not statistical acceptance thresholds.
policy={'numeric_rtol':1e-5,'numeric_atol':1e-8,'reason':'Allow small floating-point/optimizer platform differences; report every exceedance for inspection, not silently relax. Discrete labels and counts must match exactly. Large optimizer differences can change scientific conclusions and require review.'}
(out/'comparison-policy.json').write_text(json.dumps(policy,indent=2))
state={'run_directory':str(run),'state':'running','scripts':[],'started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
def save(): (out/'run-status.json').write_text(json.dumps(state,indent=2))
save()
for name in ['run.py','followup.py','circulation.py']:
 start=time.time();state['current_script']=name;save();print('START',name,flush=True)
 with (out/(name+'.stdout.log')).open('w',encoding='utf-8') as so,(out/(name+'.stderr.log')).open('w',encoding='utf-8') as se:
  q=subprocess.run([sys.executable,'-X','utf8',str(run/'companion_causal_test'/name)],cwd=run,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','PYTHONUTF8':'1'},stdout=so,stderr=se,timeout=900)
 state['scripts'].append({'script':name,'exit_code':q.returncode,'seconds':time.time()-start});save();print('FINISH',name,q.returncode,flush=True)
 if q.returncode:state['state']='script_failed';save();sys.exit(q.returncode)
comparisons=[]
def walk(a,b,path,entry):
 if isinstance(a,(int,float)) and not isinstance(a,bool) and isinstance(b,(int,float)) and not isinstance(b,bool):
  entry['numeric_leaves']+=1
  err=abs(a-b);entry['max_absolute_difference']=max(entry['max_absolute_difference'],err)
  if not math.isclose(a,b,rel_tol=policy['numeric_rtol'],abs_tol=policy['numeric_atol']):entry['differences'].append({'path':path,'original':a,'rerun':b,'absolute_difference':err})
 elif isinstance(a,dict) and isinstance(b,dict):
  if set(a)!=set(b):entry['differences'].append({'path':path,'keys_original':sorted(a),'keys_rerun':sorted(b)})
  for k in a.keys()&b.keys():walk(a[k],b[k],path+'/'+str(k),entry)
 elif isinstance(a,list) and isinstance(b,list):
  if len(a)!=len(b):entry['differences'].append({'path':path,'length_original':len(a),'length_rerun':len(b)})
  for i,(x,y) in enumerate(zip(a,b)):walk(x,y,path+'/'+str(i),entry)
 elif a!=b:entry['differences'].append({'path':path,'original':a,'rerun':b})
for name in ['results.json','followup_results.json','circulation_results.json','protocol.json']:
 a=json.loads((root/'companion_causal_test'/name).read_text(encoding='utf-8'));b=json.loads((run/'companion_causal_test'/name).read_text(encoding='utf-8'))
 entry={'file':name,'numeric_leaves':0,'max_absolute_difference':0.,'differences':[]};walk(a,b,'',entry);comparisons.append(entry)
for name in ['energy.csv','cv_predictions.csv','cluster_predictions.csv','control_predictions.csv']:
 def read(p):
  rows=list(csv.DictReader(p.open(encoding='utf-8')))
  for row in rows:
   for k,v in row.items():
    try:row[k]=float(v)
    except (ValueError,TypeError):pass
  return rows
 entry={'file':name,'numeric_leaves':0,'max_absolute_difference':0.,'differences':[]};walk(read(root/'companion_causal_test'/name),read(run/'companion_causal_test'/name),'',entry);comparisons.append(entry)
(out/'result-differences.json').write_text(json.dumps(comparisons,indent=2))
preserved=all((root/p).exists() and sha(root/p)==h for p,h in original_hashes.items())
state.update(state='finished',original_outputs_preserved=preserved,total_exceedances=sum(len(c['differences']) for c in comparisons),finished_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()));save()
print(json.dumps(state),flush=True)
