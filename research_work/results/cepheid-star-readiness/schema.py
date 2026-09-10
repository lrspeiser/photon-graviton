"""Read official archive schema before constructing the catalog query."""
from pathlib import Path
import requests

HERE=Path(__file__).resolve().parent
query="SELECT column_name, datatype, unit FROM TAP_SCHEMA.columns WHERE table_name='gaiadr3.vari_cepheid'"
response=requests.get('https://gea.esac.esa.int/tap-server/tap/sync',params={
    'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':query},timeout=45)
response.raise_for_status()
(HERE/'cepheid-schema.txt').write_text(response.text,encoding='utf-8',newline='\n')
print(response.text)
