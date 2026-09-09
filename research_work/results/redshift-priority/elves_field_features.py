"""Export only allowed identity/method fields; never parse distance/velocity values."""
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

HERE=Path(__file__).resolve().parent
EXPOSED={'IC1613','UGC00685','UGC00695','UGC01056','UGC01085','dw0020p0837','PiscesA','dw0112p0129'}


def parse(source):
    raw=source.read_bytes()
    # Publisher header ends after the note block, at its last dashed separator.
    text=raw.decode('utf-8-sig')
    marker='-'*80
    data=text.rsplit(marker,1)[1]
    records=[]
    for line in data.splitlines():
        if not line.strip():continue
        assert len(line)>=153, 'Unexpected short data row'
        name=line[:13].strip();method=line[31:32];flags=line[154:159].strip()
        # Presence only: no float conversion, export or comparison of the velocity.
        has_velocity=bool(line[144:153].strip())
        reasons=[]
        if method not in ('T','S'):reasons.append('redshift-derived distance' if method=='R' else 'unrecognized distance method')
        if not has_velocity:reasons.append('no listed velocity')
        if method=='S':
            parts=flags.split('-')
            if len(parts)!=3 or any(x not in ('T','F') for x in parts):reasons.append('unrecognized SBF flags')
            elif parts[1]!='F' or parts[2]!='F':reasons.append('ambiguous or failed SBF measurement')
        if name in EXPOSED:reasons.append('example row already exposed by web result')
        records.append(dict(target_name=name,ra_deg=float(line[14:22]),dec_deg=float(line[23:30]),
            distance_method=method,velocity_present=has_velocity,sbf_flags=flags,
            provisional_eligible=not reasons,exclusion_reasons=reasons,
            freshness='unresolved; alias, group, historical exposure and frame audit required'))
    assert len(records)==95 and len({r['target_name'] for r in records})==95
    assert EXPOSED.issubset({r['target_name'] for r in records})
    return raw,records


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('source',type=Path);args=parser.parse_args()
    raw,records=parse(args.source)
    assert hashlib.sha256(raw).hexdigest()=='890b1d178feb505cb2d9105926d13b77bc639d882e18776a3b3561b0650d8903', 'Source version changed; audit the new release before use'
    result=dict(source_url='https://content.cld.iop.org/journals/0004-637X/1001/2/244/revision1/apjae4c5ct1_mrt.txt',
        source_sha256=hashlib.sha256(raw).hexdigest(),protocol_sha256=hashlib.sha256((HERE/'elves-field-audit-protocol.md').read_bytes()).hexdigest(),
        rows=len(records),method_counts=dict(Counter(r['distance_method'] for r in records)),
        provisional_eligible=sum(r['provisional_eligible'] for r in records),
        velocity_present_by_method=dict(Counter(r['distance_method'] for r in records if r['velocity_present'])),
        explicitly_exposed_names=sorted(EXPOSED),
        status='Feature-only candidate staging; not a certified fresh sample, no target values scored.',records=records)
    (HERE/'elves-field-feature-audit.json').write_text(json.dumps(result,indent=2)+'\n',newline='\n')
    print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2))
