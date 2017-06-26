from dipsim import multiframe, util, detector, illuminator, microscope 
import numpy as np
import matplotlib.pyplot as plt

import os; import time; start = time.time(); print('Running...')

# Main input parameters
n_pts = 250 # 50000
pols = [np.array([-1,1,0]), np.array([1,1,0]), np.array([0,1,0]),np.array([1,0,0])]
ill_axes = [np.array([0, 0, 1]), np.array([0, 0, 1]), np.array([0, 0, 1]), np.array([0, 0, 1])]
det_axes = [np.array([0, 0, 1]), np.array([0, 0, 1]), np.array([-1.0/np.sqrt(2), 0, 1.0/np.sqrt(2)]), np.array([-1.0/np.sqrt(2), 0, 1.0/np.sqrt(2)])]
n_rows = 4
n_cols = len(det_axes)
inch_fig = 5
vis_px = 2000
dpi = 250
col_labels = 4*['']
row_labels = ['Scene', r'Excitation Efficiency $\eta_{\text{exc}}$', r'Detection Efficiency $\eta_{\text{det}}$', r'Total Efficiency $\eta_{\text{tot}}$']

# Generate axes
size = (inch_fig*n_cols, inch_fig*n_rows)
fig, axs = plt.subplots(n_rows, n_cols, figsize=size)
plt.subplots_adjust(wspace=0.2, hspace=-0.3)
if len(pols) == 1:
    axs = np.expand_dims(axs, 1)
caxs = util.generate_caxs(axs)


# Compute and plot on axes
for i, (pol, det_axis, ill_axis) in enumerate(zip(pols, det_axes, ill_axes)):
    print('Computing microscope: ' + str(i))

    # Create microscope
    ill = illuminator.Illuminator(illum_type='kohler',
                                  optical_axis=ill_axis,
                                  na=0.8,
                                  bfp_pol_dir=pol)
    
    det = detector.Detector(optical_axis=np.array(det_axis),
                            na=0.8,
                            n=1.5)
                                  
    m = microscope.Microscope(illuminator=ill, detector=det, max_photons=1)
    
    # Plot scene and efficiencies
    m.plot_excitation_efficiency(n=n_pts, my_ax=axs[1,i], my_cax=caxs[1,i],
                                 color_min=0, color_max=1.0)
    m.plot_collection_efficiency(n=n_pts, my_ax=axs[2,i], my_cax=caxs[2,i],
                                 color_min=0, color_max=0.06)
    m.plot_sensitivity(n=n_pts, my_ax=axs[3,i], my_cax=caxs[3,i],
                                 color_min=0, color_max=0.06)
    
    m.draw_scene(my_ax=axs[0,i], dpi=dpi, vis_px=vis_px)
    caxs[0,i].axis('off')

# Label axes and save
util.label_rows_and_cols(axs, row_labels, col_labels)
print('Saving final figure.')    
fig.savefig('detection.png', dpi=dpi)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
