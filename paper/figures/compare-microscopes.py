from dipsim import util
from dipsim import multiframe
import numpy as np
import matplotlib.pyplot as plt
import os; import time; start = time.time(); print('Running...')

# Main input parameters
row_labels = ['Geometry', r'$\sigma_{\Omega}$']
col_labels = ['Single-view (NA${}_\\textrm{ill}$=0, NA${}_\\textrm{det}$=1.1', 'Dual-view oblique symmetric widefield (NA=0.6, $\beta$=53${}^{\circ}$', 'Dual-view orthogonal symmetric light-sheet (NA=0.8)', 'Dual-view orthogonal symmetric light-sheet (NA=0.94)']

n_pts = 10000
n_cols = len(col_labels)
n_rows = 2
inch_fig = 5
dpi = 200

m1 = multiframe.MultiFrameMicroscope(ill_thetas=[0], det_thetas=[0],
                                    ill_nas=[0.8], det_nas=[1.1],
                                    ill_types=['sheet'], det_types=['lens'],
                                    colors=['(1,0,0)'], n_frames=4,
                                    n_pts=n_pts, max_photons=1000, n_samp=1.33)

angle = np.pi/2 - np.arcsin(0.6/1.33)
m2 = multiframe.MultiFrameMicroscope(ill_thetas=[angle, -angle], det_thetas=[-angle, angle],
                                    ill_nas=2*[0.6], det_nas=2*[0.6],
                                    ill_types=2*['wide'], det_types=2*['lens'],
                                    colors=['(1,0,0)','(0,0,1)'], n_frames=4,
                                    n_pts=n_pts, max_photons=500, n_samp=1.33)

m3 = multiframe.MultiFrameMicroscope(ill_thetas=[np.pi/4, -np.pi/4], det_thetas=[-np.pi/4, np.pi/4],
                                    ill_nas=2*[0], det_nas=2*[0.8],
                                    ill_types=2*['sheet'], det_types=2*['lens'],
                                    colors=['(1,0,0)','(0,0,1)'], n_frames=4,
                                    n_pts=n_pts, max_photons=500, n_samp=1.33)


m4 = multiframe.MultiFrameMicroscope(ill_thetas=[np.pi/4, -np.pi/4], det_thetas=[-np.pi/4, np.pi/4],
                                    ill_nas=2*[0], det_nas=2*[1.33*np.sin(np.pi/4)],
                                    ill_types=2*['sheet'], det_types=2*['lens'],
                                    colors=['(1,0,0)','(0,0,1)'], n_frames=4,
                                    n_pts=n_pts, max_photons=500, n_samp=1.33)

experiments = [m1, m2, m3, m4]

# Generate axes
size = (inch_fig*n_cols, inch_fig*n_rows)
fig, axs = plt.subplots(n_rows, n_cols, figsize=size)
plt.subplots_adjust(wspace=0.2, hspace=0)
if len(col_labels) == 1:
    axs = np.expand_dims(axs, 1)
caxs = util.generate_caxs(axs)

# Compute and plot on axes
for i, exp in enumerate(experiments):
    exp.calc_estimation_stats()
    
    scene_string = exp.scene_string()
    util.draw_scene(scene_string, my_ax=axs[0,i], dpi=dpi)
    
    util.plot_sphere(directions=exp.directions, data=exp.sa_uncert,
                     color_norm='log', linthresh=1e-4,
                     color_min=8e-4, color_max=1e0,
                     my_ax=axs[1,i], my_cax=caxs[1,i])

    caxs[0,i].axis('off')
    
    data = exp.sa_uncert
    maximum = np.max(data)
    med = np.median(data)
    mad = np.median(np.abs(data - med))

    print(str(col_labels[i]) + '   {:.2e}'.format(maximum) + '   {:.2e}'.format(med) + '   {:.2e}'.format(mad))

# Label axes and save
util.label_rows_and_cols(axs, row_labels, col_labels, row_pos=(-0.2, 0.5))
print('Saving final figure.')    
fig.savefig('compare-microscopes.pdf', dpi=dpi)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
