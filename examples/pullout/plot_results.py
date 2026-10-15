from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p = Path(__file__).resolve().parent
fig, axs = plt.subplots(2, 1, figsize=(8, 8), layout='constrained')
for job, colour, label in [('cdp', '#246bc1', 'CDP, damage cap 0.99'),
                            ('cdpm2', '#d52228', 'CDPM2')]:
    data = np.atleast_2d(np.loadtxt(p / (job + '_curve.csv'), delimiter=',', skiprows=1))
    if data.shape[1] != 5 or not np.isfinite(data).all():
        raise ValueError('Invalid extracted curve: ' + job)
    for ax in axs:
        ax.plot(data[:, 3], data[:, 4], color=colour, linewidth=2, label=label)
        ax.set_xlabel('Unloaded-end relative slip (mm)')
        ax.set_ylabel('Full-specimen load (kN)')
        ax.grid(alpha=0.2)
axs[0].set_title('Load-slip comparison')
axs[0].legend()
axs[1].set_title('Initial region')
axs[1].set_xlim(0, 0.6)
fig.savefig(p / 'comparison.png', dpi=180)
