from dipsim import multiframe, util
import numpy as np
import matplotlib.pyplot as plt

import os; import time; start = time.time(); print('Running...')

# Main input parameters
n_pts = 10000
bfp_n = 256
illum_det_angles = np.deg2rad([0, 45, 90])
n_cols = len(illum_det_angles)
n_rows = 3
inch_fig = 5
vis_px = 2000
dpi = 250
row_labels = ['Scene', r'$\sigma_{\Omega} = \sigma_{\Phi}\sigma_{\Theta}\sin\Theta$', r'Fractional Histogram']
col_labels = ['Epi-detection', r'$45^{\circ}$-detection', r'Ortho-detection']

# Generate axes
size = (inch_fig*n_cols, inch_fig*n_rows)
fig, axs = plt.subplots(n_rows, n_cols, figsize=size)
plt.subplots_adjust(wspace=0.2, hspace=0)
if len(illum_det_angles) == 1:
    axs = np.expand_dims(axs, 1)
caxs = util.generate_caxs(axs)

# Compute and plot on axes
for i, (illum_det_angle) in enumerate(illum_det_angles):
    print('Computing microscope: ' + str(illum_det_angle))
    m = multiframe.TwoViewPolScope(n_pts=n_pts, n_frames=4, 
                                   illum_det_angle=illum_det_angle,
                                   #na_ill=0.8, na_det=0.8, n_samp=1.33,
                                   na1=0.8, na2=0.8, n_samp=1.33,
                                   det_type='lens', bfp_n=bfp_n,
                                   dist_type='poisson', max_photons=500)

    m.draw_scene(my_ax=axs[0,i], dpi=dpi, vis_px=vis_px)
    
    util.plot_sphere(directions=m.directions, data=m.sa_uncert,
                     color_norm='log', linthresh=1e-4,
                     color_min=1e-3, color_max=4*np.pi,
                     my_ax=axs[1,i], my_cax=caxs[1,i])

    util.plot_histogram(m.sa_uncert, ax=axs[2,i])

    caxs[0,i].axis('off')
    caxs[2,i].axis('off')        

# Label axes and save
util.label_rows_and_cols(axs, row_labels, col_labels, row_pos=(-0.2, 0.5))
print('Saving final figure.')    
fig.savefig('compare_hist.png', dpi=dpi)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
