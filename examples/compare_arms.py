from dipsim import multiframe, util
import numpy as np
import matplotlib.pyplot as plt

import os; import time; start = time.time(); print('Running...')

# Main input parameters
n_pts = 100000
bfp_n = 256
illum_det_angles = np.deg2rad([0, 45, 90])
gauss_stds = [0]#,1,3]
n_cols = len(illum_det_angles)
n_rows = 5
inch_fig = 5
vis_px = 2000
dpi = 250
row_labels = ['Scene', r'$F^{-1}_{0,0} = \mathrm{Var}(\Theta) = \sigma_{\Theta}^2$', r'$F^{-1}_{1,1} = \mathrm{Var}(\Phi) = \sigma_{\Phi}$', r'$F^{-1}_{0,1} = \mathrm{Cov}(\Theta, \Phi)$', r'$\sigma_{\Omega} = \sigma_{\Phi}\sigma_{\Theta}\sin\Theta$']
col_labels = ['Epi-detection', r'$45^{\circ}$-detection', r'Ortho-detection']

# Generate axes
size = (inch_fig*n_cols, inch_fig*n_rows)
fig, axs = plt.subplots(n_rows, n_cols, figsize=size)
plt.subplots_adjust(wspace=0.2, hspace=-0.3)
if len(illum_det_angles) == 1:
    axs = np.expand_dims(axs, 1)
caxs = util.generate_caxs(axs)

# Compute and plot on axes
for i, (illum_det_angle) in enumerate(illum_det_angles):
    print('Computing microscope: ' + str(illum_det_angle))
    m = multiframe.OneArmPolScope(n=n_pts, n_frames=4, 
                                  illum_det_angle=illum_det_angle,
                                  det_type='lens', bfp_n=bfp_n,
                                  dist_type='poisson', max_photons=1000)

    data = [m.fi_inv[:,0], m.fi_inv[:,3], m.fi_inv[:,1], m.sa_uncert]
    scales = ['log', 'log', 'linlog', 'log']
    for j, (d, scale) in enumerate(zip(data, scales)):
        util.plot_sphere(directions=m.directions, data=d,
                         color_norm=scale, linthresh=1e-4,
                         color_min=1e-3, color_max=4*np.pi,
                         my_ax=axs[j+1,i], my_cax=caxs[j+1,i])
                     
    m.draw_scene(my_ax=axs[0,i], dpi=dpi, vis_px=vis_px)
    caxs[0,i].axis('off')    

# Label axes and save
util.label_rows_and_cols(axs, row_labels, col_labels)
print('Saving final figure.')    
fig.savefig('compare_arms.png', dpi=dpi)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
