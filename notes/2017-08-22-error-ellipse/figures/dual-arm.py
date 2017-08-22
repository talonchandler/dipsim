from dipsim import multiframe, util
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as patches
import os; import time; start = time.time(); print('Running...')
import matplotlib.gridspec as gridspec

# Main input parameters
col_labels = ['Geometry\n(NA${}_{\\textrm{1}}$ = 0.8, NA${}_{\\textrm{2}}$ = 0.8)', 'Uncertainty Ellipses', r'$\sigma_{\Omega}$ [sr]']
fig_labels = ['a)', 'b)', 'c)']
n_pts = 100 # Points on sphere
n_pts_sphere = 50000 # Points on sphere
n_grid_pts = 5
inch_fig = 5
dpi = 300

# Setup figure and axes
fig = plt.figure(figsize=(3.2*inch_fig, 1*inch_fig))
gs0 = gridspec.GridSpec(1, 3, wspace=0.2, hspace=0.1)
gs00 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[0,0], width_ratios=[1, 0.05], wspace=0.1)
gs10 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[0,1], width_ratios=[1, 0.05], wspace=0.1)
gs20 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[0,2], width_ratios=[1, 0.05], wspace=0.1)

ax0 = plt.subplot(gs00[0])
ax1 = plt.subplot(gs10[0])
ax2 = plt.subplot(gs20[0])
cax2 = plt.subplot(gs20[1])

for ax, col_label, fig_label  in zip([ax0, ax1, ax2], col_labels, fig_labels):
    ax.annotate(col_label, xy=(0,0), xytext=(0.5, 1.05), textcoords='axes fraction',
                va='bottom', ha='center', fontsize=14, annotation_clip=False)
    ax.annotate(fig_label, xy=(0,0), xytext=(0, 1.05), textcoords='axes fraction',
                va='bottom', ha='center', fontsize=14, annotation_clip=False)
    
for ax in [ax0, ax1, ax2]:
    ax.tick_params(axis='both', labelsize=14)
for cax in [cax2]:
    cax.tick_params(axis='both', labelsize=14)

# Calculate a list of points to sample in region
n = 1.33
NA_ill = np.linspace(0, 1.33, num=1000)
NA_det = NA_ill
pts = np.mgrid[n/n_grid_pts/2:n-n/n_grid_pts/2:n_grid_pts*1j,n/n_grid_pts/2:n-n/n_grid_pts/2:n_grid_pts*1j].reshape((2, n_grid_pts**2)).T.tolist()

def is_feasible(pt):
    if pt[0] < pt[1] + 0.2:
        return True
    else:
        return False

pts_list = [pt for pt in pts if is_feasible(pt)]
pts = np.array(pts_list).T

exp = multiframe.MultiFrameMicroscope(ill_thetas=[np.pi/4, -np.pi/4], det_thetas=[-np.pi/4, np.pi/4],
                                      ill_nas=[0.8, 0.8], det_nas=[0.8, 0.8],
                                      ill_types=2*['sheet'], det_types=2*['lens'],
                                      colors=['(1,0,0)', '(0,0,1)'], n_frames=4,
                                      n_pts=n_pts_sphere, max_photons=[1, 1], n_samp=1.33)

exp.calc_estimation_stats()

util.draw_scene(exp.scene_string(), my_ax=ax0, dpi=dpi)
util.draw_scene(exp.ellipse_string(n_pts=250), my_ax=ax1, dpi=dpi)
util.plot_sphere(directions=exp.directions, data=exp.sa_uncert,
                 color_norm='log', linthresh=1e-4,
                 color_min=1e-4, color_max=1e1,
                 my_ax=ax2, my_cax=cax2)
    
# Label axes and save
print('Saving final figure.')    
fig.savefig('dual-view.pdf', dpi=1000)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
