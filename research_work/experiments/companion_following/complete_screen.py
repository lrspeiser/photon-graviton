"""Complete only the unavailable-in-batch ring screen; preserve first execution."""
import json
from pathlib import Path
from model import catalogue, Mechanics, initial
from campaign import HERE, manifest, simulate, write_json

def main():
    evidence = HERE/"evidence"
    original = evidence/"screen-v1"
    output = evidence/"screen-completion-v1"
    output.mkdir(exist_ok=False)
    write_json(output/"manifest.json", manifest())
    good=[]; unavailable=[]
    for row in catalogue()[600:]:
        try:
            initial(Mechanics([row]), "ring")
            good.append(row)
        except ValueError as exc:
            unavailable.append(dict(id=row["id"], fixture="ring", duration=64,
                behavioral_pass=False, conservation_pass=False, circulation_min=0.,
                unavailable=True, reason=str(exc), trajectory_executed=False))
    write_json(output/"unavailable.json", unavailable)
    rings = simulate(good, "ring", output/"conservative-ring")
    screen = {
        "phenomenological": {f: json.loads((original/f"phenomenological-{f}"/"summary.json").read_text())["results"]
                            for f in ("chain","ring")},
        "conservative": {"chain": json.loads((original/"conservative-chain"/"summary.json").read_text())["results"],
                         "ring": rings+unavailable}}
    write_json(output/"screen.json", screen)
    write_json(output/"complete.json", dict(complete=True, unchanged_first_execution=str(original),
               unavailable_cases=len(unavailable), ring_trajectories=len(good)))
    print(f"Completed original screen without changing a fixture: {len(good)} rings; {len(unavailable)} unavailable.")

if __name__=="__main__": main()
