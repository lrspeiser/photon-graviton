from pathlib import Path
import re,json
P=Path(__file__).resolve().parent
canon=lambda s:re.sub(r'(?<=\D)0+(?=\d)','',re.sub(r'[^A-Z0-9]','',s.upper()))
dust={}
for l in (P/'data/themis.dat').read_text().splitlines():
 try:
  d=dict(name=l[:23].strip(),ra=float(l[24:33]),dec=float(l[34:43]),dist=float(l[49:61]),lum=float(l[100:109]),elum=float(l[110:119]),mass=float(l[80:89]),sfr=float(l[62:70]))
  if d['lum']>0 and d['dist']>0:dust[canon(d['name'])]=d
 except ValueError:pass
split=json.loads((P.parent/'companion_wave_test/data/sparc_frozen.json').read_text())['split']
for s,ns in split.items():print(s,len(ns),sum(canon(n) in dust for n in ns))
print('Dust rows',len(dust));print([(n,dust[canon(n)]['name']) for ns in split.values() for n in ns if canon(n) in dust])
