from dipsim import multiframe, util, detector, illuminator, microscope, util
import numpy as np
import matplotlib.pyplot as plt

import os; import time; start = time.time(); print('Running...')

# Main input parameters
n_pts = 250
pols = 4*[np.pi/4]
ill_axes = 4*[np.pi/4]
det_axes = 4*[-np.pi/4]
ill_nas = 4*[0.8]
det_nas = 4*[0.8]
kappas = [-3, 0, 3, 6]
n_rows = len(det_axes)
n_cols = 5
inch_fig = 5
dpi = 300

col_labels = ['Fluorophore Geometry', 'Microscope Geometry', r'Absorption Efficiency $\eta_{\textrm{abs}}$', r'Detection Efficiency $\eta_{\textrm{det}}$', r'Total Efficiency $\eta_{\textrm{tot}}$']
row_labels = ['$\kappa = -3$', '$\kappa = 0$', '$\kappa = 3$', '$\kappa = 6$']

# Generate axes
size = (inch_fig*n_cols, inch_fig*n_rows)
fig, axs = plt.subplots(n_rows, n_cols, figsize=size)
plt.subplots_adjust(wspace=0.1, hspace=-0.1)
if len(pols) == 1:
    axs = np.expand_dims(axs, 1)
caxs = util.generate_caxs(axs)

# Compute and plot on axes
for i, (pol, det_axis, ill_axis, det_na, ill_na, kappa) in enumerate(zip(pols, det_axes, ill_axes, det_nas, ill_nas, kappas)):
    print('Computing microscope: ' + str(i))

    # Create microscope
    ill = illuminator.Illuminator(illum_type='wide',
                                  theta_optical_axis=ill_axis,
                                  na=ill_na,
                                  phi_pol=pol)
    
    det = detector.Detector(theta_optical_axis=np.array(det_axis),
                            na=det_na)
                                  
    m = microscope.Microscope(illuminator=ill, detector=det, max_photons=1)
    
    # Plot scene and efficiencies
    m.plot_excitation_efficiency(n=n_pts, kappa=kappa, my_ax=axs[i,2], my_cax=caxs[i,2],
                                 color_min=0, color_max=0.75)
    m.plot_collection_efficiency(n=n_pts, kappa=kappa, my_ax=axs[i,3], my_cax=caxs[i,3],
                                 color_min=0, color_max=0.25)
    m.plot_sensitivity(n=n_pts, kappa=kappa, my_ax=axs[i,4], my_cax=caxs[i,4],
                                 color_min=0, color_max=0.15)

    scene_string = m.scene_string()
    util.draw_scene(scene_string, my_ax=axs[i,1], dpi=dpi)
    caxs[i,1].axis('off')

    fluorophore_string = 'watson(0, 0, kappa, 0, 0, 0);\n'
    fluorophore_string = fluorophore_string.replace('kappa', str(kappa))
    util.draw_scene(fluorophore_string, my_ax=axs[i,0], dpi=dpi, chop=False)
    caxs[i,0].axis('off')

# Label axes and save
util.label_rows_and_cols(axs, row_labels, col_labels)
print('Saving final figure.')    
fig.savefig('single-frame.pdf', dpi=dpi)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
