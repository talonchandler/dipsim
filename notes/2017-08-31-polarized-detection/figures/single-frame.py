from dipsim import multiframe, util, detector, illuminator, microscope, util
import numpy as np
import matplotlib.pyplot as plt

import os; import time; start = time.time(); print('Running...')

# Main input parameters
n_pts = 1000
ill_types = ['unpolarized', 'unpolarized', 'wide']
ill_pols = [0, np.pi/4, np.pi/2]
det_pols = [0, np.pi/4, np.pi/2]
ill_axes = [0, -np.pi/4, np.pi/4]
det_axes = [0, np.pi/4, np.pi]
ill_nas = [0.8, 0.8, 0.8]
det_nas = [1.0, 0.8, 0.8]
n_rows = len(det_axes)
n_cols = 4
inch_fig = 5
dpi = 300

col_labels = ['Geometry', r'Absorption Efficiency $\eta_{\textrm{abs}}$', r'Detection Efficiency $\eta_{\textrm{det}}$', r'Total Efficiency $\eta_{\textrm{tot}}$']
row_labels = 3*['']

# Generate axes
size = (inch_fig*n_cols, inch_fig*n_rows)
fig, axs = plt.subplots(n_rows, n_cols, figsize=size)
plt.subplots_adjust(wspace=0.1, hspace=-0.1)
if len(ill_pols) == 1:
    axs = np.expand_dims(axs, 1)
caxs = util.generate_caxs(axs)

# Compute and plot on axes
for i, (ill_type, det_pol, ill_pol, det_axis, ill_axis, det_na, ill_na) in enumerate(zip(ill_types, det_pols, ill_pols, det_axes, ill_axes, det_nas, ill_nas)):
    print('Computing microscope: ' + str(i))

    # Create microscope
    ill = illuminator.Illuminator(illum_type=ill_type,
                                  theta_optical_axis=ill_axis,
                                  na=ill_na,
                                  phi_pol=ill_pol)
    
    det = detector.Detector(det_type='polarized',
                            theta_optical_axis=np.array(det_axis),
                            na=det_na,
                            phi_pol=det_pol)
                                  
    m = microscope.Microscope(illuminator=ill, detector=det, max_photons=1)
    
    # Plot scene and efficiencies
    m.plot_excitation_efficiency(n=n_pts, my_ax=axs[i,1], my_cax=caxs[i,1],
                                 color_min=0, color_max=1.0)
    m.plot_collection_efficiency(n=n_pts, my_ax=axs[i,2], my_cax=caxs[i,2],
                                 color_min=0, color_max=0.25)
    m.plot_sensitivity(n=n_pts, my_ax=axs[i,3], my_cax=caxs[i,3],
                                 color_min=0, color_max=0.25)

    scene_string = m.scene_string()
    util.draw_scene(scene_string, my_ax=axs[i,0], dpi=dpi)
    caxs[i,0].axis('off')

# Label axes and save
util.label_rows_and_cols(axs, row_labels, col_labels)
print('Saving final figure.')    
fig.savefig('single-frame.pdf', dpi=dpi)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
