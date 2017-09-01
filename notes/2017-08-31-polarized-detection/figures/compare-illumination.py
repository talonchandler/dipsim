from dipsim import util
from dipsim import multiframe
import numpy as np
import matplotlib.pyplot as plt
import os; import time; start = time.time(); print('Running...')

# Main input parameters
col_labels = ['Geometry', 'Uncertainty Ellipse', r'$\sigma_{\Omega}$']
row_labels = 6*['']

n_pts = 10000
n_rows = len(row_labels)
n_cols = 3
inch_fig = 5
dpi = 200

m1 = multiframe.MultiFrameMicroscope(ill_thetas=[0], det_thetas=[0],
                                    ill_nas=[0.01], det_nas=[0.8],
                                    ill_types=['sheet'], det_types=['lens'],
                                    colors=['(1,0,0)'], n_frames=2,
                                    n_pts=n_pts, max_photons=4000/2, n_samp=1.33)

m2 = multiframe.MultiFrameMicroscope(ill_thetas=[0], det_thetas=[0],
                                    ill_nas=[0.01], det_nas=[0.8],
                                    ill_types=['sheet'], det_types=['lens'],
                                    colors=['(1,0,0)'], n_frames=3,
                                    n_pts=n_pts, max_photons=4000/3, n_samp=1.33)

m3 = multiframe.MultiFrameMicroscope(ill_thetas=[0], det_thetas=[0],
                                    ill_nas=[0.01], det_nas=[0.8],
                                    ill_types=['sheet'], det_types=['lens'],
                                    colors=['(1,0,0)'], n_frames=4,
                                    n_pts=n_pts, max_photons=4000/4, n_samp=1.33)

m4 = multiframe.MultiFrameMicroscope(ill_thetas=[np.pi/4, -np.pi/4], det_thetas=[-np.pi/4, np.pi/4],
                                     ill_nas=2*[0.01], det_nas=2*[0.8],
                                     ill_types=2*['sheet'], det_types=2*['lens'],
                                     colors=['(1,0,0)','(0,0,1)'], n_frames=2,
                                     n_pts=n_pts, max_photons=4000/4, n_samp=1.33)

m5 = multiframe.MultiFrameMicroscope(ill_thetas=[np.pi/4, -np.pi/4], det_thetas=[-np.pi/4, np.pi/4],
                                     ill_nas=2*[0.01], det_nas=2*[0.8],
                                     ill_types=2*['sheet'], det_types=2*['lens'],
                                     colors=['(1,0,0)','(0,0,1)'], n_frames=3,
                                     n_pts=n_pts, max_photons=4000/6, n_samp=1.33)

m6 = multiframe.MultiFrameMicroscope(ill_thetas=[np.pi/4, -np.pi/4], det_thetas=[-np.pi/4, np.pi/4],
                                     ill_nas=2*[0.01], det_nas=2*[0.8],
                                     ill_types=2*['sheet'], det_types=2*['lens'],
                                     colors=['(1,0,0)','(0,0,1)'], n_frames=4,
                                     n_pts=n_pts, max_photons=4000/8, n_samp=1.33)

experiments = [m1, m2, m3, m4, m5, m6]

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
    
    util.draw_scene(exp.scene_string(), my_ax=axs[i,0], dpi=dpi)
    util.draw_scene(exp.ellipse_string(n_pts=250), my_ax=axs[i,1], dpi=dpi)
    
    util.plot_sphere(directions=exp.directions, data=exp.sa_uncert,
                     color_norm='log', linthresh=1e-4,
                     color_min=1e-4, color_max=1e1,
                     my_ax=axs[i,2], my_cax=caxs[i,2])

    caxs[i,0].axis('off')
    caxs[i,1].axis('off')
    
    data = exp.sa_uncert
    maximum = np.max(data)
    med = np.median(data)
    mad = np.median(np.abs(data - med))

    print(str(row_labels[i]) + '   {:.2e}'.format(maximum) + '   {:.2e}'.format(med) + '   {:.2e}'.format(mad))

# Label axes and save
util.label_rows_and_cols(axs, row_labels, col_labels, row_pos=(-0.2, 0.5))
print('Saving final figure.')    
fig.savefig('compare-illumination.pdf', dpi=dpi)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
