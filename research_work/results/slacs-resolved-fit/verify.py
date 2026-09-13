from pathlib import Path
import json,importlib.util
import numpy as np
from model import AnnularModel
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('orbit',HERE.parent/'slacs-orbit-sensitivity/model.py');orbit=importlib.util.module_from_spec(spec);spec.loader.exec_module(orbit)
model=AnnularModel(3.,np.array([0.,1.,2.,4.]),1.5,.24,.46,22000.,20.)
checks=[]
for beta in [-.3,0,.3]:
 ann=model.coefficients(beta)
 combined=(ann*model.den).sum(axis=1)/model.den.sum()
 aperture=orbit.coefficients(3.,4.,1.5,.24,.46,22000.,20.,[beta])[0]
 error=float(np.max(abs(combined/aperture-1)));assert error<1e-10
 checks.append(dict(beta=beta,relative_partition_error=error))
(HERE/'verification.json').write_text(json.dumps(dict(scope='Flux-weighted annular second moments reconstruct the previously checked cumulative aperture',checks=checks),indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(checks))
