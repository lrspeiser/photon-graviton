"""Additional local thin-source audit; does not alter production quadrature."""
import json
import numpy as np
from run import HERE, Model

results=[]
for R in (30,100):
    for nq in (8,16,32):
        model=Model(R,256,nq)
        k=1e-10
        shell,_,_=model.absorption(np.full(model.n,k))
        error=float(max(abs(shell/(k*model.V)-1)))
        results.append(dict(R=R,impact_nodes=nq,thin_uniform_local_source_max_relative_error=error,
                            passes_supplementary_gate=error<.001))
(HERE/'thin-source-check.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf8',newline='\n')
assert all(x['passes_supplementary_gate'] for x in results if x['impact_nodes']>=16)
