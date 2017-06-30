from dipsim import util
from dipsim import multiframe
import numpy as np
import matplotlib.pyplot as plt
import os; import time; start = time.time(); print('Running...')

# Main input parameters
row_labels = ['Geometry', r'$\sigma_{\Omega} = \sqrt{\textrm{det}\{I^{-1}\}}\sin\Theta$', r'Fractional Histogram']
col_labels = ['One-View Wide-Field \n0.8 NA', 'Orthogonal Two-View Wide-Field \n0.8/0.8 NA', 'Two-View Wide-Field \n0.8/0.8 NA', 'Symmetric Two-View Light-Sheet \n0.8/0.8 NA', 'Asymmetric Two-View Light-Sheet \n0.71/1.1 NA', 'Three-View Wide-Field \n0.8/0.8/1.2 NA', 'Three-View Light-Sheet \n0.8/0.8/1.2 NA']

n_pts = 50000
n_cols = len(col_labels)
n_rows = 3
inch_fig = 5
dpi = 400

m1 = multiframe.MultiFrameMicroscope(ill_thetas=[0], det_thetas=[0],
                                    ill_nas=[0.8], det_nas=[0.8],
                                    ill_types=['wide'], det_types=['lens'],
                                    colors=['(1,0,0)'], n_frames=4,
                                    n_pts=n_pts, max_photons=1000, n_samp=1.33)

m2 = multiframe.MultiFrameMicroscope(ill_thetas=[np.pi/4, -np.pi/4], det_thetas=[-np.pi/4, np.pi/4],
                                    ill_nas=2*[0.8], det_nas=2*[0.8],
                                    ill_types=2*['wide'], det_types=2*['lens'],
                                    colors=['(1,0,0)','(0,0,1)'], n_frames=4,
                                    n_pts=n_pts, max_photons=500, n_samp=1.33)

m3 = multiframe.MultiFrameMicroscope(ill_thetas=[np.pi/4, -np.pi/4], det_thetas=[np.pi/4, -np.pi/4],
                                    ill_nas=2*[0.8], det_nas=2*[0.8],
                                    ill_types=2*['wide'], det_types=2*['lens'],
                                    colors=['(1,0,0)','(0,0,1)'], n_frames=4,
                                    n_pts=n_pts, max_photons=500, n_samp=1.33)

m4 = multiframe.MultiFrameMicroscope(ill_thetas=[np.pi/4, -np.pi/4], det_thetas=[-np.pi/4, np.pi/4],
                                    ill_nas=2*[0], det_nas=2*[0.8],
                                    ill_types=2*['sheet'], det_types=2*['lens'],
                                    colors=['(1,0,0)','(0,0,1)'], n_frames=4,
                                    n_pts=n_pts, max_photons=500, n_samp=1.33)

m5 = multiframe.MultiFrameMicroscope(ill_thetas=[np.deg2rad(57), np.deg2rad(-33)], det_thetas=[np.deg2rad(-33), np.deg2rad(57)],
                                    ill_nas=2*[0], det_nas=[1.1, 0.71],
                                    ill_types=2*['sheet'], det_types=2*['lens'],
                                    colors=['(1,0,0)','(0,0,1)'], n_frames=4,
                                    n_pts=n_pts, max_photons=500, n_samp=1.33)

m6 = multiframe.MultiFrameMicroscope(ill_thetas=[np.pi/4, -np.pi/4, np.pi/4, -np.pi/4], det_thetas=[-np.pi/4, np.pi/4, np.pi, np.pi],
                                    ill_nas=4*[0.8], det_nas=[0.8, 0.8, 1.2, 1.2],
                                    ill_types=4*['wide'], det_types=4*['lens'],
                                    colors=['(1,0,0)','(0,0,1)','(1,0,0)','(0,0,1)'], n_frames=4,
                                    n_pts=n_pts, max_photons=500, n_samp=1.33)


m7 = multiframe.MultiFrameMicroscope(ill_thetas=[np.pi/4, -np.pi/4, np.pi/4, -np.pi/4], det_thetas=[-np.pi/4, np.pi/4, np.pi, np.pi],
                                    ill_nas=4*[0], det_nas=[0.8, 0.8, 1.2, 1.2],
                                    ill_types=4*['sheet'], det_types=4*['lens'],
                                    colors=['(1,0,0)','(0,0,1)','(1,0,0)','(0,0,1)'], n_frames=4,
                                    n_pts=n_pts, max_photons=500, n_samp=1.33)

experiments = [m1, m2, m3, m4, m5, m6, m7]

# Generate axes
size = (inch_fig*n_cols, inch_fig*n_rows)
fig, axs = plt.subplots(n_rows, n_cols, figsize=size)
plt.subplots_adjust(wspace=0.2, hspace=0)
if len(col_labels) == 1:
    axs = np.expand_dims(axs, 1)
caxs = util.generate_caxs(axs)

# Compute and plot on axes
for i, exp in enumerate(experiments):
    print('Computing microscope: ' + str(col_labels[i]))
    exp.calc_estimation_stats()
    
    scene_string = exp.scene_string()
    util.draw_scene(scene_string, my_ax=axs[0,i], dpi=dpi)
    
    util.plot_sphere(directions=exp.directions, data=exp.root_det_sin,
                     color_norm='log', linthresh=1e-4,
                     color_min=8e-4, color_max=1e0,
                     my_ax=axs[1,i], my_cax=caxs[1,i])

    util.plot_histogram(exp.root_det_sin, ax=axs[2,i])

    caxs[0,i].axis('off')
    caxs[2,i].axis('off')        

# Label axes and save
util.label_rows_and_cols(axs, row_labels, col_labels, row_pos=(-0.2, 0.5))
print('Saving final figure.')    
fig.savefig('compare-microscopes.png', dpi=dpi)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
