from dipsim import multiframe, util
import numpy as np
import matplotlib.pyplot as plt

import os; import time; start = time.time(); print('Running...')

# Main input parameters
n_pts = 5000 # 50000
bfp_n = 64
n_frames = [3, 4]
noise_type = ['poisson', 'gaussian', 'gaussian', 'gaussian']
gauss_stds = [0, 0.5, 1, 3]
n_cols = len(gauss_stds) + 1
n_rows = len(n_frames)
inch_fig = 5
vis_px = 2000
dpi = 500
col_labels = ['Scene', 'Poisson', 'P+Gaussian: $\eta=0, \sigma=0.5$', 'P+Gaussian: $\eta=0, \sigma=1.0$', 'P+Gaussian: $\eta=0, \sigma=3.0$']
row_labels = ['N=3', 'N=4']

# Generate axes
size = (inch_fig*n_cols, inch_fig*n_rows)
fig, axs = plt.subplots(n_rows, n_cols, figsize=size)
plt.subplots_adjust(wspace=0.2, hspace=0)
if n_cols == 1:
    axs = np.expand_dims(axs, 1)
caxs = util.generate_caxs(axs)

# Compute and plot on axes
for j, nf in enumerate(n_frames):
    for i, std in enumerate(gauss_stds):
        print('Computing microscope: ' + str(std))    
        m = multiframe.NFramePolScope(n_frames=nf, bfp_n=bfp_n, dist_type=noise_type[i], gauss_mean=0, gauss_std=std)
        m.plot_orientation_std(filename=str(i)+'frame.png', n=n_pts, color_norm='log',
                               my_ax=axs[j,i+1], my_cax=caxs[j,i+1],
                               color_min=1e-3, color_max=4*np.pi)
        m.draw_scene(my_ax=axs[j,0], dpi=dpi, vis_px=vis_px)
        caxs[j,0].axis('off')

# Label axes and save
util.label_rows_and_cols(axs, row_labels, col_labels)
print('Saving final figure.')    
fig.savefig('compare_noise_models.png', dpi=dpi)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
