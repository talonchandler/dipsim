from dipsim import multiframe, util
import numpy as np
import matplotlib.pyplot as plt

import os; import time; start = time.time(); print('Running...')

# Main input parameters
n_pts = 100000
bfp_n = 256
det_angles = np.deg2rad([0, 0, 30, 60, 85, 89, 90])
#det_angles = det_angles + [x*np.pi/6 for x in range(5)]
#det_axes = [np.array([0,0,1]), np.array([0,0,1]), np.array([-1,0,0]), np.array([-1.0/np.sqrt(2),0,1.0/np.sqrt(2)])]
det_types = ['4pi', 'lens', 'lens', 'lens', 'lens', 'lens', 'lens']
noise_types = ['poisson']#, 'gaussian', 'gaussian']
gauss_stds = [0]#,1,3]
n_cols = len(det_angles)
n_rows = len(noise_types) + 1
inch_fig = 5
vis_px = 2000
dpi = 250
row_labels = ['Scene', r'Poisson']#, r'Poisson + Gaussian: $\sigma = 0.5$', r'Poisson + Gaussian: $\sigma = 3$']
col_labels = ['4$\pi$ detection', 'Epi-detection', r'$30^{\circ}$-detection', r'$60^{\circ}$-detection', r'$85^{\circ}$-detection', r'$89^{\circ}$-detection', r'$90^{\circ}$-detection']

# Generate axes
size = (inch_fig*n_cols, inch_fig*n_rows)
fig, axs = plt.subplots(n_rows, n_cols, figsize=size)
plt.subplots_adjust(wspace=0.2, hspace=0)
if len(det_angles) == 1:
    axs = np.expand_dims(axs, 1)
caxs = util.generate_caxs(axs)

# Compute and plot on axes
for i, (det_angle, det_type) in enumerate(zip(det_angles, det_types)):
    for j, (noise_type, gauss_std) in enumerate(zip(noise_types, gauss_stds)):    
        print('Computing microscope: ' + str(det_angle))    
        m = multiframe.NFramePolScope(n_frames=4, det_axis=np.array([-np.sin(det_angle), 0, np.cos(det_angle)]),
                                      det_type=det_type, bfp_n=bfp_n,
                                      dist_type=noise_type, gauss_mean=0,
                                      gauss_std=gauss_std,
                                      max_photons=1000)                                      
                                      # crlb_frame=util.rot_mat(np.pi/4, np.array([0,1,0])),

        m.plot_orientation_std(filename=str(i)+'frame.png', n=n_pts, color_norm='log',
                               my_ax=axs[j+1,i], my_cax=caxs[j+1,i],
                               color_min=1e-3, color_max=4*np.pi)
    m.draw_scene(my_ax=axs[0,i], dpi=dpi, vis_px=vis_px)
    caxs[0,i].axis('off')    

# Label axes and save
util.label_rows_and_cols(axs, row_labels, col_labels)
print('Saving final figure.')    
fig.savefig('compare_detection.png', dpi=dpi)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
