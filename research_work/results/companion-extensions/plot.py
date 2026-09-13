from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ref=json.loads((P/'refraction-results.json').read_text());phase=json.loads((P/'phase-results.json').read_text())
fig,axes=plt.subplots(2,2,figsize=(12,8),sharex=True)
for j,b in enumerate(['I','II']):
 runs=[q for q in ref['runs'] if q['baryons']==b]
 control=next(q for q in runs if q['epsilon_infinity']==1)
 ph=next(q for q in phase['selected'] if q['baryons']==b)
 curves=[('Exact-third control',control,'#333333'),('Phase selector',ph,'#008c95')]
 for driver,color in [('ordinary','#e68923'),('companion','#9557b4')]:
  best=min((q for q in runs if q['driver']==driver),key=lambda q:q['scores']['inner']['RMSE_kms'])
  curves.append((driver.title()+'-driven refraction',best,color))
 axes[0,j].scatter(Robs,Vobs,s=16,c='black',label='Published circular-speed bins',zorder=5)
 for label,q,color in curves:
  v=np.array(q['predicted_kms']);axes[0,j].plot(Robs,v,label=label,c=color);axes[1,j].plot(Robs,v-Vobs,c=color)
 for ax in axes[:,j]:
  ax.axvspan((Robs[19]+Robs[20])/2,26,color='#c5d8ef',alpha=.3);ax.grid(alpha=.2);ax.set_xlim(5,26)
 axes[0,j].set_title('Ordinary-matter baseline '+b);axes[1,j].axhline(0,c='black',lw=.6);axes[1,j].set_xlabel('Galactocentric radius (kpc)')
axes[0,0].set_ylabel('Circular speed (km/s)');axes[1,0].set_ylabel('Prediction minus observed (km/s)');axes[0,0].legend(fontsize=8)
fig.suptitle('Fixed one-third capture with proposed extensions\nShaded bins excluded from parameter selection; previously seen data. Gas-monopole pilot.',fontsize=12)
fig.tight_layout();fig.savefig(P/'comparison.png',dpi=160)
