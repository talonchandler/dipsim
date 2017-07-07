from dipsim import util
from dipsim import multiframe
import numpy as np
import matplotlib.pyplot as plt
import os; import time; start = time.time(); print('Running...')

# Main input parameters
row_labels = ['Geometry', r'$\sigma_{\Omega}$', r'Relative Frequency']
nas = [0.7, 0.9, 1.1, 1.3]
col_labels = ['NA='+str(na) for na in nas]

n_pts = 500
n_cols = len(col_labels)
n_rows = 3
inch_fig = 5
dpi = 400

# Generate axes
size = (inch_fig*n_cols, inch_fig*n_rows)
fig, axs = plt.subplots(n_rows, n_cols, figsize=size)
plt.subplots_adjust(wspace=0.2, hspace=0)
if len(col_labels) == 1:
    axs = np.expand_dims(axs, 1)
caxs = util.generate_caxs(axs)

# Compute and plot on axes
for i, na in enumerate(nas):
    exp = multiframe.MultiFrameMicroscope(ill_thetas=[0], det_thetas=[0],
                                          ill_nas=[na], det_nas=[na],
                                          ill_types=['wide'], det_types=['lens'],
                                          colors=['(1,0,0)'], n_frames=4,
                                          n_pts=n_pts, max_photons=1000, n_samp=1.33)
    
    print('Computing microscope: ' + str(col_labels[i]))
    exp.calc_estimation_stats()
    
    scene_string = exp.scene_string()
    util.draw_scene(scene_string, my_ax=axs[0,i], dpi=dpi)
    
    util.plot_sphere(directions=exp.directions, data=exp.sa_uncert,
                     color_norm='log', linthresh=1e-4,
                     color_min=8e-4, color_max=1e0,
                     my_ax=axs[1,i], my_cax=caxs[1,i])

    util.plot_histogram(exp.sa_uncert, ax=axs[2,i])

    caxs[0,i].axis('off')
    caxs[2,i].axis('off')        

# Label axes and save
util.label_rows_and_cols(axs, row_labels, col_labels, row_pos=(-0.2, 0.5))
print('Saving final figure.')    
fig.savefig('single-crlb.png', dpi=dpi)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
