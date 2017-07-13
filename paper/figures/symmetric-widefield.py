from dipsim import multiframe, util
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as patches
import os; import time; start = time.time(); print('Running...')
import matplotlib.gridspec as gridspec

# Main input parameters
col_labels = ['Geometry', r'$\sigma_{\Omega}$', 'Mean$\{\sigma_{\Omega}\}$', 'STD$\{\sigma_{\Omega}\}$']
n_pts = 100 # Points on sphere
n_pts_sphere = 10000 # Points on sphere
n_grid_pts = 13
n_rows, n_cols = 1, len(col_labels)
inch_fig = 5
dpi = 300

# Setup figure and axes
fig = plt.figure(figsize=(4.5*inch_fig, inch_fig))
gs0 = gridspec.GridSpec(1, 2, wspace=0.125)
gs00 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs0[0], width_ratios=[1, 1, 0.05], wspace=0.05)
gs01 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[1], width_ratios=[1, 1], wspace=0.3)
gs010 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs01[0], width_ratios=[1, 0.05], wspace=0.1)
gs011 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs01[1], width_ratios=[1, 0.05], wspace=0.1)
ax0 = plt.subplot(gs00[0])
ax1 = plt.subplot(gs00[1])
ax2 = plt.subplot(gs010[0])
ax3 = plt.subplot(gs011[0])
cax1 = plt.subplot(gs00[2])
cax2 = plt.subplot(gs010[1])
cax3 = plt.subplot(gs011[1])

for ax, col_label in zip([ax0, ax1, ax2, ax3], col_labels):
    ax.annotate(col_label, xy=(0,0), xytext=(0.5, 1.1), textcoords='axes fraction',
                va='center', ha='center', fontsize=18, annotation_clip=False)

# Calculate a list of points to sample in region
n = 1.33
NA_max = n*np.sin(np.pi/4)
NA = np.linspace(0, n, 1000)
lens_bound = np.rad2deg(2*np.arcsin(NA/n))
cover_bound = np.rad2deg(np.pi - 2*np.arcsin(NA/n))

pts = np.mgrid[n/n_grid_pts/2:n:n_grid_pts*1j,0:180:n_grid_pts*1j].reshape((2, n_grid_pts**2)).T.tolist()

def is_feasible(pt):
    if pt[1] < np.rad2deg(np.pi - 2*np.arcsin(pt[0]/n)) + 8 and pt[1] > np.rad2deg(2*np.arcsin(pt[0]/n)) - 8 and pt[0] < NA_max + 0.1:
        return True
    else:
        return False

pts_list = [pt for pt in pts if is_feasible(pt)]
pts = np.array(pts_list).T

# Calculate mean and variance for each point
def calc_stats(param):
    na = param[0]
    angle = param[1]
    exp = multiframe.MultiFrameMicroscope(ill_thetas=[np.deg2rad(angle/2)], det_thetas=[-np.deg2rad(angle/2)],
                                      ill_nas=[na], det_nas=[na],
                                      ill_types=['wide'], det_types=['lens'],
                                      colors=['(1,0,0)'], n_frames=4,
                                      n_pts=n_pts, max_photons=500, n_samp=1.33)

    exp.calc_estimation_stats()
    return np.mean(exp.sa_uncert), np.std(exp.sa_uncert)

mean = []
std = []
for i, pt in enumerate(pts.T):
    print('Calculating microscope '+str(i+1)+'/'+str(pts.shape[1]))
    x = calc_stats(pt)
    mean.append(x[0])
    std.append(x[1])

# Plot 2D regions
def plot_2d_regions(ax, cax, pts, data):
    ax.plot(NA, lens_bound, 'k-', zorder=11)
    ax.plot(NA, cover_bound, 'k-', zorder=11)

    # Set y ticks
    from matplotlib.ticker import FuncFormatter, FixedLocator
    def degrees(x, pos):
        return str(int(x)) + '${}^{\circ}$'
    ax.yaxis.set_major_locator(FixedLocator([0, 45, 90, 135, 180]))
    ax.yaxis.set_major_formatter(FuncFormatter(degrees))

    # Annotation
    def my_annotate(ax, annotation, xy, fontsize=9, rotation=0):
        ax.annotate(annotation, xy=(0,0), xytext=xy, textcoords='axes fraction',
                    va='center', ha='center', fontsize=fontsize,
                    annotation_clip=False, rotation=rotation, zorder=13)

    my_annotate(ax, 'NA', (0.5, -0.10), fontsize=12)
    my_annotate(ax, 'Angle Between Arms', (-0.125, 0.5), fontsize=12, rotation=90)
    my_annotate(ax, 'Objectives collide\nwith cover slip', (0.65, 0.85), fontsize=12)
    my_annotate(ax, 'Objectives collide\nwith each other', (0.65, 0.15), fontsize=12)
    my_annotate(ax, 'Feasible', (0.3, 0.5), fontsize=12)

    # Calculate colors
    color_map='coolwarm'
    color_norm='log'
    color_min=np.min(data)#1e-3
    color_max=np.max(data)#1e2
    if color_norm == 'linear':
        norm = matplotlib.colors.Normalize(vmin=color_min, vmax=color_max)
    elif color_norm == 'log':
        norm = matplotlib.colors.LogNorm(vmin=color_min, vmax=color_max)
    elif color_norm == 'linlog':
        norm = matplotlib.colors.SymLogNorm(linthresh=linthresh, vmin=-color_max, vmax=color_max)
    elif color_norm == 'power':
        norm = matplotlib.colors.PowerNorm(gamma=gamma, vmin=data.min(), vmax=data.max())

    norm_data = norm(data).data
    norm_data2 = np.expand_dims(norm_data, 1)    
    cmap = matplotlib.cm.get_cmap(color_map)
    colors = np.apply_along_axis(cmap, 1, norm_data2)

    # Plot scatter for colorbar
    sc = ax.scatter(pts[0,:], pts[1,:], c=data, s=0, cmap=cmap, norm=norm,
                    marker='s', lw=0)

    # Plot patches
    width = n/(n_grid_pts)
    for i, (pt, c) in enumerate(zip(pts_list, colors)):
        if pt[1] == 0:
            height = 180/14.5
        if pt[1] == 0:
            height = 180/(n_grid_pts-0.5)
        ax.add_patch(patches.Rectangle((pt[0] - width/2, pt[1] - height/2), width, height, facecolor=c, edgecolor=c))
    fig.colorbar(sc, cax=cax, orientation='vertical')

    # Mask around lines
    ax.fill_between(NA, lens_bound, 0, color='white', zorder=10)
    ax.fill_between(NA, cover_bound, 180, color='white', zorder=10)

    ax.set(xlim=[0, 1.33], ylim=[0, 180])
    
# Plot first two columns
angle = 60
na = 0.6
exp = multiframe.MultiFrameMicroscope(ill_thetas=[np.deg2rad(angle/2)], det_thetas=[-np.deg2rad(angle/2)],
                                      ill_nas=[na], det_nas=[na],
                                      ill_types=['wide'], det_types=['lens'],
                                      colors=['(1,0,0)'], n_frames=4,
                                      n_pts=n_pts_sphere, max_photons=500, n_samp=1.33)
exp.calc_estimation_stats()
scene_string = exp.scene_string()

util.draw_scene(scene_string, my_ax=ax0, dpi=dpi)
util.plot_sphere(directions=exp.directions, data=exp.sa_uncert,
                 color_norm='log', linthresh=1e-4,
                 color_min=8e-4, color_max=1e0,
                 my_ax=ax1, my_cax=cax1)
    
# Plots last two columns
plot_2d_regions(ax2, cax2, pts, mean)
plot_2d_regions(ax3, cax3, pts, std)
    
# Label axes and save
axs = np.array([ax0, ax1, ax2, ax3])
print('Saving final figure.')    
fig.savefig('symmetric-widefield.pdf', dpi=250)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')