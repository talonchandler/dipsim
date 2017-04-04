from dipsim import multiframe, util
import numpy as np
import matplotlib.pyplot as plt

import os; import time; start = time.time(); print('Running...')

# Main input parameters
n_pts = 50000
n_frames = [1,1] #n_frames = [1, 2, 3, 4]
n_cols = len(n_frames) 
n_rows = 7
inch_fig = 5
vis_px = 2000
dpi = 500
row_labels = ['Scene', r'$\sigma_{\phi}^*$', r'$\sigma_{\theta}^*$',r'$\sigma_{\Omega} = \sigma_{\theta}^*\sigma_{\phi}^*$', r'$\frac{dI}{d\theta^*}$', r'$\frac{dI}{d\phi^*}$', r'$I$']
col_labels = [r'$N='+str(x)+'$' for x in n_frames]

# Generate axes
size = (inch_fig*n_cols, inch_fig*n_rows)
fig, axs = plt.subplots(n_rows, n_cols, figsize=size)
plt.subplots_adjust(wspace=0.2, hspace=0)
if len(n_frames) == 1:
    axs = np.expand_dims(axs, 1)
caxs = util.generate_caxs(axs)

# Compute and plot on axes
for i, nf in enumerate(n_frames):
    print('Computing microscope: ' + str(nf))    
    m = multiframe.NFramePolScope(n_frames=nf)
    m.plot_orientation_std(n=n_pts, color_norm='linlog',
                           my_axs=axs[1:,i], my_caxs=caxs[1:,i])
    m.draw_scene(my_ax=axs[0,i], dpi=dpi, vis_px=vis_px)
    caxs[0,i].axis('off')

# Label axes and save
util.label_rows_and_cols(axs, row_labels, col_labels)
print('Saving final figure.')    
fig.savefig('compare_N_pol.png', dpi=dpi)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
