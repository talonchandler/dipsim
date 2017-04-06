from dipsim import multiframe, util
import numpy as np
import matplotlib.pyplot as plt

import os; import time; start = time.time(); print('Running...')

# Main input parameters
n_pts = 50 # 50000
bfp_n = 4 # 256
types = ['poisson', 'gaussian', 'gaussian', 'gaussian']
stds = [0, 0.5, 1, 3]
n_cols = 1
n_rows = 1
inch_fig = 5
vis_px = 2000
dpi = 500
row_labels = ['Detector Term']

# Generate axes
size = (inch_fig*n_cols, inch_fig*n_rows)
fig, axs = plt.subplots(n_rows, n_cols, figsize=size)
plt.subplots_adjust(wspace=0.2, hspace=0.4)
axs = np.expand_dims(axs, 1)
caxs = util.generate_caxs(axs)

# Compute and plot on axes
for i, std in enumerate(stds):
    print('Computing microscope: ' + str(i))    
    m = multiframe.NFramePolScope(n_frames=1, bfp_n=bfp_n, dist_type=types[i], gauss_mean=0, gauss_std=std)
    m.noise_model.plot_detector_term(my_ax=axs[0])
    # m.noise_model.plot_pdf(my_ax=axs[1])    
    
    caxs[0].axis('off')
    # caxs[1].axis('off')    

# Label axes and save
print('Saving final figure.')    
fig.savefig('noise_model_test.png', dpi=dpi)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
