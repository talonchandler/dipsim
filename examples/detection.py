from dipsim import multiframe, util, detector, illuminator, microscope 
import numpy as np
import matplotlib.pyplot as plt

import os; import time; start = time.time(); print('Running...')

# Main input parameters
n_pts = 5000 # 50000
bfp_n = 64
n_frames = [4]
n_cols = 3
n_rows = len(n_frames)
inch_fig = 5
vis_px = 2000
dpi = 500
row_labels = ['N=4']
col_labels = ['Scene', 'Excitation Efficiency', 'Collection Efficiency']

# Generate axes
size = (inch_fig*n_cols, inch_fig*n_rows)
fig, axs = plt.subplots(n_rows, n_cols, figsize=size)
plt.subplots_adjust(wspace=0.2, hspace=0.4)
if len(n_frames) == 1:
    axs = np.expand_dims(axs, 1)
caxs = util.generate_caxs(axs)


# Compute and plot on axes
for i, nf in enumerate(n_frames):
    print('Computing microscope: ' + str(nf))
    # Constant detection path.
    det = detector.Detector(optical_axis=np.array([0, 0, 1]),
                            na=1.3,
                            n=1.5)
    # Changing illumination.
    
    ill = illuminator.Illuminator(illum_type='kohler',
                                  optical_axis=np.array([0., 0., 1.]),
                                  f=10,
                                  bfp_rad=2,
                                  bfp_pol=np.array([1,0,0]),
                                  bfp_n=bfp_n)
    m = microscope.Microscope(illuminator=ill, detector=det, max_photons=1)
    m.plot_excitation_efficiency(n=n_pts, my_ax=axs[1,i], my_cax=caxs[1,i],
                                 color_min=0, color_max=1.0)
    m.plot_collection_efficiency(n=n_pts, my_ax=axs[2,i], my_cax=caxs[2,i],
                                 color_norm='linear',
                                 color_min=0.05, color_max=0.1)

    m.draw_scene(my_ax=axs[0,i], dpi=dpi, vis_px=vis_px)
    caxs[0,i].axis('off')

# Label axes and save
util.label_rows_and_cols(axs, row_labels, col_labels)
print('Saving final figure.')    
fig.savefig('detection.png', dpi=dpi)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
