"""Correct lag-window scoring from saved trajectories; preserve every first score."""
import copy
import json
from model import leader_angle
from campaign import HERE, manifest, write_json, best_correlation
from archive_io import open_npz
import numpy as np

def main():
    evidence=HERE/"evidence"
    root=evidence/"scoring-v2"
    root.mkdir(exist_ok=False)
    write_json(root/"manifest.json",manifest())
    source_screen=json.loads((evidence/"screen-completion-v1"/"screen.json").read_text())
    screen=copy.deepcopy(source_screen)
    scores={};changed=[]
    for summary in sorted(evidence.rglob("summary.json")):
        data=json.loads(summary.read_text())
        if data["fixture"] != "chain": continue
        raw=open_npz(summary.with_name("trajectories.npz"))
        t=raw["time"]
        lead=np.array([leader_angle(tt,data["amplitude"]) for tt in t])
        for i,r in enumerate(data["results"]):
            state=raw["state"][i]
            signal=np.mean(np.arctan2(state[:,-4:,3],state[:,-4:,2]),axis=1)
            corr,lag=best_correlation(signal,lead)
            before=dict(pass_value=r["behavioral_pass"],correlation=r["last_four_correlation"],lag=r["best_lag"])
            r["last_four_correlation"]=corr;r["best_lag"]=lag
            r["behavioral_pass"]=bool(data["amplitude"] != 0 and r["chain_gain"]>=.1
                and corr>=.7 and r["maximum_follower_gain"]<=3
                and r["min_separation"]>.1 and r["first_nonfinite_time"] is None)
            if before["pass_value"]!=r["behavioral_pass"] or abs(before["correlation"]-corr)>1e-12:
                changed.append(dict(source=str(summary.relative_to(evidence)),id=r["id"],before=before,
                     after=dict(pass_value=r["behavioral_pass"],correlation=corr,lag=lag)))
        scores[summary.parent.relative_to(evidence).as_posix()]=data
        for group in ("phenomenological","conservative"):
            if summary.parent==evidence/"screen-v1"/f"{group}-chain":
                screen[group]["chain"]=data["results"]
    write_json(root/"screen.json",screen)
    write_json(root/"rescored-prior-chains.json",scores)
    write_json(root/"changes.json",changed)
    write_json(root/"complete.json",dict(complete=True,trajectories_reintegrated=0))
    for g in screen:
        a=screen[g]["chain"];b={r["id"]:r for r in screen[g]["ring"]}
        print(g,"corrected chain passes",sum(r["behavioral_pass"] for r in a),
              "joint passes",sum(r["behavioral_pass"] and b[r["id"]]["behavioral_pass"] for r in a))
if __name__=="__main__": main()
