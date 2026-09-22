"""Restricted FITS IMAGE reader, no external FITS dependency.

Supports only primary/IMAGE HDUs with PCOUNT=0,GCOUNT=1 and standard numeric
BITPIX values. Rejects other layouts, including compressed/tabled images.
The original .fits.gz bytes are gzip streams containing uncompressed HDUs.
No arbitrary FITS expression evaluation is used.
"""
import gzip
from pathlib import Path
import numpy as np


def _value(raw):
    raw=raw.strip()
    if raw.startswith("'"):
        i=1;s=''
        while i<len(raw):
            if raw[i]=="'":
                if i+1<len(raw) and raw[i+1]=="'":s+="'";i+=2;continue
                return s.strip()
            s+=raw[i];i+=1
        raise ValueError('Unterminated FITS string')
    raw=raw.split('/')[0].strip()
    if raw in ('T','F'):return raw=='T'
    try:return int(raw)
    except ValueError:
        try:return float(raw.replace('D','E'))
        except ValueError:return raw


def read_images(path):
    path=Path(path);data=gzip.open(path,'rb').read() if path.suffix=='.gz' else path.read_bytes()
    pos=0;out={};idx=0
    while pos<len(data):
        if len(data)-pos<2880 or not data[pos:pos+80].strip(b' \0'):break
        start=pos;header={};found=False
        while pos+80<=len(data):
            card=data[pos:pos+80].decode('ascii');pos+=80;k=card[:8].strip()
            if k=='END':found=True;break
            if card[8:10]=='= ':header[k]=_value(card[10:])
        if not found:raise ValueError('No FITS END')
        pos=start+((pos-start+2879)//2880)*2880
        if header.get('XTENSION','IMAGE')!='IMAGE':raise ValueError('Not a supported IMAGE HDU')
        if header.get('PCOUNT',0)!=0 or header.get('GCOUNT',1)!=1:raise ValueError('Unsupported group layout')
        bit=int(header['BITPIX']);dtype={8:'u1',16:'>i2',32:'>i4',64:'>i8',-32:'>f4',-64:'>f8'}[bit]
        shape=tuple(int(header[f'NAXIS{j}']) for j in range(int(header['NAXIS']),0,-1))
        count=int(np.prod(shape)) if shape else 0;size=count*(abs(bit)//8)
        if pos+size>len(data):raise ValueError('Truncated FITS image')
        a=np.frombuffer(data,dtype=dtype,count=count,offset=pos).reshape(shape).copy() if count else None
        if a is not None and (header.get('BSCALE',1)!=1 or header.get('BZERO',0)!=0):
            a=a.astype(float)*header.get('BSCALE',1)+header.get('BZERO',0)
        name=header.get('EXTNAME','PRIMARY' if idx==0 else f'HDU{idx}')
        out[name]=(header,a);pos+=((size+2879)//2880)*2880;idx+=1
    return out


def channel(header,name):
    found=[int(k[1:])-1 for k,v in header.items() if k.startswith('C') and k[1:].isdigit() and str(v).strip()==name]
    if len(found)!=1:raise ValueError(f'Channel {name} not unique: {found}')
    return found[0]
