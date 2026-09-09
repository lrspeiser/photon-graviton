from pathlib import Path
import urllib.request,urllib.parse,json
P=Path(__file__).resolve().parent
params={'format':'json','COMMAND':"'-82'",'EPHEM_TYPE':"'VECTORS'",'CENTER':"'500@399'",'START_TIME':"'2005-01-01'",'STOP_TIME':"'2006-01-01'",'STEP_SIZE':"'30 d'",'OUT_UNITS':"'KM-S'",'VEC_TABLE':"'3'",'VEC_CORR':"'NONE'",'CSV_FORMAT':"'YES'"}
url='https://ssd.jpl.nasa.gov/api/horizons.api?'+urllib.parse.urlencode(params)
result=json.loads(urllib.request.urlopen(url,timeout=40).read())
(P/'data/horizons_response.json').write_text(json.dumps(result,indent=2));(P/'data/horizons_request.json').write_text(json.dumps({'url':url,'parameters':params},indent=2))
print(result.get('result',result)[:1600])
