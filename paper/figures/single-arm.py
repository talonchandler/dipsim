from dipsim import util
from dipsim import multiframe
import numpy as np
import matplotlib.pyplot as plt
import os; import time; start = time.time(); print('Running...')

# Main input parameters
row_labels = ['']
col_labels = ['Geometry', r'$\sigma_{\Omega}$']

n_pts = 100000
n_cols = len(col_labels)
n_rows = 1
inch_fig = 5
dpi = 300

# Generate axes
size = (inch_fig*n_cols, inch_fig*n_rows)
fig, axs = plt.subplots(n_rows, n_cols, figsize=size)
plt.subplots_adjust(wspace=0, hspace=0)
if len(col_labels) == 1:
    axs = np.expand_dims(axs, 1)
if len(row_labels) == 1:
    axs = np.expand_dims(axs, 0)
caxs = util.generate_caxs(axs)

# Compute and plot on axes
exp = multiframe.MultiFrameMicroscope(ill_thetas=[0], det_thetas=[0],
                                      ill_nas=[0.8], det_nas=[0.8],
                                      ill_types=['wide'], det_types=['lens'],
                                      colors=['(1,0,0)'], n_frames=4,
                                      n_pts=n_pts, max_photons=1000, n_samp=1.33)

print('Computing microscope')
exp.calc_estimation_stats()

scene_string = exp.scene_string()
util.draw_scene(scene_string, my_ax=axs[0,0], dpi=dpi)

util.plot_sphere(directions=exp.directions, data=exp.sa_uncert,
                 color_norm='log', linthresh=1e-4,
                 color_min=8e-4, color_max=1e0,
                 my_ax=axs[0,1], my_cax=caxs[0,1])

caxs[0,0].axis('off')

# Label axes and save
util.label_rows_and_cols(axs, row_labels, col_labels, row_pos=(-0.2, 0.5))
print('Saving final figure.')    
fig.savefig('../paper/single-arm.pdf', dpi=dpi)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
