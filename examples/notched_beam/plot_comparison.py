# SPDX-License-Identifier: LGPL-2.1-or-later
"""Plot the six locally extracted analysis curves."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
here=Path(__file__).resolve().parent
curves=here/'results'
fig,axes=plt.subplots(1,3,figsize=(11,3.5),sharey=True)
for ax,size,count in zip(axes,['coarse','medium','fine'],[70,280,1120]):
    for model,color in [('cdp','#2472ad'),('cdpm2','#d34a35')]:
        a=np.loadtxt(curves/(model+'_'+size+'_curve.txt'))
        ax.plot(a[:,0],a[:,1],color=color,lw=2,label=model.upper())
    ax.set(title=str(count)+' elements',xlabel='Displacement (mm)',xlim=(0,1.55),ylim=(0,.95))
    ax.grid(alpha=.2);ax.legend(fontsize=6)
axes[0].set_ylabel('Full-beam load at 50 mm width (kN)')
fig.tight_layout();folder=here/'results';folder.mkdir(exist_ok=True)
fig.savefig(folder/'comparison.png',dpi=180)
plt.close(fig)
