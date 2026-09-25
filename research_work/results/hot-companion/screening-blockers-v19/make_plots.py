#!/usr/bin/env python3
"""Plot the already-computed model checks. No astronomical observations are plotted."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)

def rows(name):
    with (ROOT/'results'/name).open() as f:return list(csv.DictReader(f))

rr=rows('stochastic_aggregate.csv')
r=[d for d in rr if d['stage']=='launched' and float(d['time_tau'])==1.]
x=np.array([float(d['x']) for d in r]);p=np.array([float(d['simulated_open']) for d in r]);se=np.array([float(d['standard_error']) for d in r])
xx=np.linspace(0,10,500)
fig,ax=plt.subplots(figsize=(7.4,4.7))
ax.plot(xx,np.exp(-xx)*(1-np.exp(-1)),label='Derived constant-field response at t = τ')
ax.errorbar(x,p,yerr=2*se,fmt='o',capsize=3,label='Event simulation, ±2 sampling standard errors')
ax.set_yscale('log');ax.set_xlabel('Background field / screening scale, g / g_d');ax.set_ylabel('Released fraction')
ax.set_title('One million realizations per condition: stochastic model, not data')
ax.legend(fontsize=8);fig.tight_layout();fig.savefig(OUT/'constant_field_check.png',dpi=170);plt.close(fig)

ss=np.linspace(0,7,700)
fig,ax=plt.subplots(figsize=(7.4,4.7))
ax.plot(ss,np.exp(-10*np.exp(-ss)),label='Individual blockers disappear')
ax.plot(ss,1+(np.exp(-10)-1)*np.exp(-ss),linestyle='--',label='Whole gate is refreshed')
ax.set_xlabel('Time after field change / common lifetime, t / τ');ax.set_ylabel('Released fraction')
ax.set_title('Same steady screening; different response to g / g_d: 10 → 0')
ax.legend();fig.tight_layout();fig.savefig(OUT/'field_step_discriminator.png',dpi=170);plt.close(fig)

fig,ax=plt.subplots(figsize=(7.4,4.7))
xx=np.linspace(.05,5,600);x=10
y=1-np.exp(-xx)
ax.plot(xx,np.exp(-x)*y,label='Equilibrated background + one launch blocker')
ax.plot(xx,y*np.exp(-x*y),linestyle='--',label='One launch blocker; no equilibrated background')
ax.set_yscale('log');ax.set_xlabel('Age / blocker lifetime, t / τ');ax.set_ylabel('Released fraction')
ax.set_title('Initial preparation is a physical assumption, not a detail')
ax.legend(fontsize=8);fig.tight_layout();fig.savefig(OUT/'initial_preparation_control.png',dpi=170);plt.close(fig)
