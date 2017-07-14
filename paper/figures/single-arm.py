from dipsim import multiframe, util
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker as ticker
import os; import time; start = time.time(); print('Running...')
import matplotlib.gridspec as gridspec

# Main input parameters
col_labels = ['Geometry (NA = 0.8)', r'$\sigma_{\Omega}$ [sr]', '', '']
fig_labels = ['a)', 'b)', 'c)', 'd)']
n_pts = 500 # Points on sphere
n_pts_sphere = 25000 # Points on sphere
n_grid_pts = 30
inch_fig = 5
dpi = 300
percentile = 99

# Setup figure and axes
fig = plt.figure(figsize=(2.2*inch_fig, 2*inch_fig))
gs0 = gridspec.GridSpec(2, 2, wspace=0.2, hspace=0.1)
gs00 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[0,0], width_ratios=[1, 0.05], wspace=0.1)
gs10 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[1,0], width_ratios=[1, 0.05], wspace=0.1)
gs01 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[0,1], width_ratios=[1, 0.05], wspace=0.1)
gs11 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[1,1], width_ratios=[1, 0.05], wspace=0.1)
ax0 = plt.subplot(gs00[0])
ax1 = plt.subplot(gs01[0])
ax2 = plt.subplot(gs10[0])
ax3 = plt.subplot(gs11[0])
cax0 = plt.subplot(gs00[1]); cax0.axis('off');
cax1 = plt.subplot(gs01[1])
cax2 = plt.subplot(gs10[1]); cax2.axis('off');
cax3 = plt.subplot(gs11[1]); cax3.axis('off');

for ax, col_label, fig_label  in zip([ax0, ax1, ax2, ax3], col_labels, fig_labels):
    ax.annotate(col_label, xy=(0,0), xytext=(0.5, 1.05), textcoords='axes fraction',
                va='center', ha='center', fontsize=14, annotation_clip=False)
    ax.annotate(fig_label, xy=(0,0), xytext=(0, 1.05), textcoords='axes fraction',
                va='center', ha='center', fontsize=14, annotation_clip=False)

for ax in [ax0, ax1, ax2, ax3]:
    ax.tick_params(axis='both', labelsize=14)
for cax in [cax1, cax2, cax3]:
    cax.tick_params(axis='both', labelsize=14)

# Calculate a list of points to sample in region
n = 1.33
nas = np.linspace(0.01, 1.33, num=n_grid_pts)

# Calculate mean and variance for each point
mean = []
std = []
for i, na in enumerate(nas):
    print('Calculating microscope '+str(i+1)+'/'+str(len(nas)))
    exp = multiframe.MultiFrameMicroscope(ill_thetas=[0], det_thetas=[0],
                                      ill_nas=[na], det_nas=[na],
                                      ill_types=['wide'], det_types=['lens'],
                                      colors=['(1,0,0)'], n_frames=4,
                                      n_pts=n_pts, max_photons=1000, n_samp=1.33)

    exp.calc_estimation_stats()
    data = exp.sa_uncert[exp.sa_uncert < np.percentile(exp.sa_uncert, percentile)]
    
    mean.append(np.mean(data))
    std.append(np.std(data))

# Plot 1D regions
def plot_1d_regions(ax, nas, data, special_pt=(0.8,0), y_string=''):

    # Annotation
    def my_annotate(ax, annotation, xy, fontsize=9, rotation=0):
        ax.annotate(annotation, xy=(0,0), xytext=xy, textcoords='axes fraction',
                    va='center', ha='center', fontsize=fontsize,
                    annotation_clip=False, rotation=rotation, zorder=13)

    my_annotate(ax, 'NA', (0.5, -0.10), fontsize=14)
    my_annotate(ax, y_string, (-0.18, 0.5), fontsize=14, rotation=90)

    # Plot
    ax.plot(nas, data, 'k-')
    ax.plot(special_pt[0], special_pt[1], 'kx', markersize=5, zorder=14)
    ax.set_yscale('log')
    ax.set(xlim=[0, 1.33], ylim=[1e-3, 1e1])
    ax.xaxis.set_ticks([0, 0.25, 0.5, 0.75, 1.00, 1.33])
    ax.set_xticklabels(['0', '0.25', '0.5', '0.75', '1.0', '1.33'])

# Plot first two columns
na = 0.8
exp = multiframe.MultiFrameMicroscope(ill_thetas=[0], det_thetas=[0],
                                  ill_nas=[na], det_nas=[na],
                                  ill_types=['wide'], det_types=['lens'],
                                  colors=['(1,0,0)'], n_frames=4,
                                  n_pts=n_pts_sphere, max_photons=1000, n_samp=1.33)

exp.calc_estimation_stats()
scene_string = exp.scene_string()
util.draw_scene(scene_string, my_ax=ax0, dpi=dpi)
util.plot_sphere(directions=exp.directions, data=exp.sa_uncert,
                 color_norm='log', linthresh=1e-4,
                 color_min=None, color_max=None,
                 my_ax=ax1, my_cax=cax1)
    
# Plots last row
plot_1d_regions(ax2, nas, mean, special_pt=(na, 0.0054), y_string='Mean${}_{0-99}$$\{\sigma_{\Omega}\}$ [sr]')
plot_1d_regions(ax3, nas, std, special_pt=(na, 0.0092), y_string='STD${}_{0-99}$$\{\sigma_{\Omega}\}$ [sr]')

# Label axes and save
print('Saving final figure.')    
fig.savefig('../paper/single-arm.pdf', dpi=250)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
