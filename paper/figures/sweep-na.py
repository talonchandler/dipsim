from dipsim import util
from dipsim import multiframe
import numpy as np
import matplotlib.pyplot as plt
import os; import time; start = time.time(); print('Running...')

# Main input parameters
n_pts = [100, 500, 1000]
nas = np.linspace(0.01, 1.33, num=50)
inch_fig = 5
dpi = 400

# Compute and plot on axes
fig, ax = plt.subplots(1, 1, figsize=(5,5))
for n_pt in n_pts:
    d = []
    for na in nas:
        exp = multiframe.MultiFrameMicroscope(ill_thetas=[0], det_thetas=[0],
                                        ill_nas=[na], det_nas=[na],
                                        ill_types=['wide'], det_types=['lens'],
                                        colors=['(1,0,0)'], n_frames=4,
                                        n_pts=n_pt, max_photons=1000, n_samp=1.33)

        print('Computing microscope: ' + str(n_pt) + ', ' + str(na))
        exp.calc_estimation_stats()
        d.append(util.dispersion_index(exp.root_det_sin))
    print('Optimal NA is: ' + str(nas[np.argmin(d)]))
    ax.plot(nas, d, '-', label='N = '+str(n_pt))
    if n_pt == n_pts[-1]:
        x = nas[np.argmin(d)]
        y = np.min(d)
        ax.arrow(x, y+0.05, 0, -0.04)
        ax.annotate('NA* = {:.3}'.format(x), xy=(x, y+0.08), xytext=(x, y+0.08), horizontalalignment='center')

ax.set_xlabel('NA')
ax.set_ylabel('Dispersion Index')
ax.set_xlim([0, 1.33])
ax.set_ylim([0, 1])
ax.legend(loc='upper right')
fig.savefig('sweep-na.png', dpi=dpi)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
