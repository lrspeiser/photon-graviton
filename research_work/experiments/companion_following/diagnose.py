"""Post-screen diagnoses declared in diagnostic-plan.md; no selection retuning."""
import json
from model import catalogue
from campaign import HERE, manifest, simulate, write_json

def main():
    root=HERE/"evidence"/"diagnostics-v2"
    root.mkdir(exist_ok=False)
    write_json(root/"manifest.json",manifest())
    selection=json.loads((HERE/"evidence"/"refine-v2"/"selection.json").read_text())
    allrows={r["id"]:r for r in catalogue()}
    p=[allrows[i] for i in selection["phenomenological"]]
    c=[allrows[i] for i in selection["conservative"]]
    base=[allrows["F0-K0-G0-M0"]]
    runs={}
    for name,fixture,kwargs in [
        ("baseline-chain","chain",{}),
        ("baseline-ring","ring",{}),
        ("baseline-long-ring","ring",dict(duration=256,dt=.02)),
        ("baseline-long-chain","chain",dict(n=24)),
        ("baseline-long-chain-extended","chain",dict(n=24,duration=64))]:
        runs[name]=simulate(base,fixture,root/name,coupling=0.,**kwargs)
    runs["selected-long-chain-extended"]=simulate(p,"chain",root/"selected-long-chain-extended",n=24,duration=64)
    runs["conservative-long-ring-dt001"]=simulate(c,"ring",root/"conservative-long-ring-dt001",duration=256,dt=.01)
    write_json(root/"diagnostics.json",runs)
    write_json(root/"complete.json",dict(complete=True))
if __name__=="__main__": main()
