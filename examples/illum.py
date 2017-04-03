from dipsim import fluorophore, illuminator, detector, microscope, stats, util, multiframe
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import time; start = time.time(); print('Running...')

from mpl_toolkits.axes_grid1 import make_axes_locatable

n = 5000
dpi = 500
vis_px = 2000
n_frames = [1,1]
#n_frames = [1, 2, 3, 4]
n_cols = len(n_frames) 
n_rows = 7
inch_fig = 5

fig, axs = plt.subplots(n_rows, n_cols,
                        figsize=(inch_fig*n_cols, inch_fig*n_rows))

caxs = []
for ax in axs.flatten():
    divider = make_axes_locatable(ax)
    caxs.append(divider.append_axes("right", size="5%", pad=0.15))
caxs = np.array(caxs).reshape(axs.shape)
    
if len(n_frames) == 1:
    axs = np.expand_dims(axs, 1)

plt.subplots_adjust(wspace=0.2, hspace=0)

for i, nf in enumerate(n_frames):
    m = multiframe.NFramePolScope(n_frames=nf)
    m.plot_orientation_std('solid'+str(i), n=n, interact=False,
                           color_norm='linlog',
                           my_axs=axs[1:,i], my_caxs=caxs[1:,i],
                           save_file=True)
    
    m.draw_scene('scene'+str(i)+'.png', interact=False, my_ax=axs[0,i], dpi=dpi, vis_px=vis_px, save_file=False)
    caxs[0,i].axis('off')

row_labels = ['Scene', r'$\sigma_{\phi}^*$', r'$\sigma_{\theta}^*$',r'$\sigma_{\Omega} = \sigma_{\theta}^*\sigma_{\phi}^*$', r'$\frac{dI}{d\theta^*}$', r'$\frac{dI}{d\phi^*}$', r'$I$']
for i, label in enumerate(row_labels):
    axs[i][0].annotate(label, xy=(0,0), xytext=(-0.1, 0.5), textcoords='axes fraction',
                       va='center', ha='center', rotation=90, fontsize=18)

#col_labels = [r'$N=1$', r'$N=2$', r'$N=3$', r'$N=4$']
col_labels = [r'$N=1$', r'$N=1$']
for i, label in enumerate(col_labels):
    axs[0][i].annotate(label, xy=(0,0), xytext=(0.5, 1.1), textcoords='axes fraction',
                       va='center', ha='center', fontsize=18)

fig.savefig('compare_N_pol.png', dpi=dpi)
print('Total time: '+str(np.round(time.time() - start, 2)))

import os
os.system('say "done"')
