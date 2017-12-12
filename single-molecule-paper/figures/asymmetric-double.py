from dipsim import multiframe, util
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as patches
import os; import time; start = time.time(); print('Running...')
import matplotlib.gridspec as gridspec

# Main input parameters
col_labels = ['Geometry\n($\\alpha_{1}$ = 60${}^{\circ}$,\n Exposure${}_1$/Exposure${}_2$ = 5)', 'Uncertainty Ellipses', r'$\sigma_{\Omega}$ [sr]', 'Median$\{\sigma_{\Omega}\}$ [sr]', 'MAD$\{\sigma_{\Omega}\}$ [sr]']
fig_labels = ['a)', 'b)', 'c)', 'd)', 'e)']
n_pts = 500 # Points on sphere
n_pts_sphere = 50000 # Points on sphere
n_grid_pts = 25
inch_fig = 5
dpi = 300

# Setup figure and axes
fig = plt.figure(figsize=(2.6*inch_fig, 2*inch_fig))
gs0 = gridspec.GridSpec(2, 1, wspace=0, hspace=0.1, height_ratios=[1,1])
gs_up = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=gs0[0], width_ratios=[1, 1, 1, 0.06], wspace=0.2)
gs_down = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[1], width_ratios=[1, 1], wspace=0.45)
gs_down_left = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_down[0], width_ratios=[1, 0.05], wspace=0.55)
gs_down_right = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_down[1], width_ratios=[1, 0.05], wspace=0.55)
ax0 = plt.subplot(gs_up[0])
ax1 = plt.subplot(gs_up[1])
ax2 = plt.subplot(gs_up[2])
cax2 = plt.subplot(gs_up[3])
ax3 = plt.subplot(gs_down_left[0])
cax3 = plt.subplot(gs_down_left[1])
ax4 = plt.subplot(gs_down_right[0])
cax4 = plt.subplot(gs_down_right[1])

for ax, col_label, fig_label  in zip([ax0, ax1, ax2, ax3, ax4], col_labels, fig_labels):
    ax.annotate(col_label, xy=(0,0), xytext=(0.5, 1.05), textcoords='axes fraction',
                va='bottom', ha='center', fontsize=14, annotation_clip=False)
    ax.annotate(fig_label, xy=(0,0), xytext=(0, 1.05), textcoords='axes fraction',
                va='bottom', ha='center', fontsize=14, annotation_clip=False)
    
for ax in [ax0, ax1, ax2, ax3, ax4]:
    ax.tick_params(axis='both', labelsize=14)
for cax in [cax2, cax3, cax4]:
    cax.tick_params(axis='both', labelsize=14)

# Calculate a list of points to sample in region
n = 1.33
alpha_max = 90
dose_width = 4
x, step = np.linspace(-2, 2, num=n_grid_pts, retstep=True, endpoint=False)
x += step/2
y, step = np.linspace(0, alpha_max, num=n_grid_pts, retstep=True, endpoint=False)
y += step/2
pts = np.array(np.meshgrid(x, y)).reshape(2, n_grid_pts**2).T

def is_feasible(pt):
    return True

pts_list = [pt for pt in pts if is_feasible(pt)]
pts = np.array(pts_list).T

# Calculate med and mad for each point
def calc_stats(param):
    na1 = n*np.sin(np.deg2rad(param[1]))
    na2 = n*np.sin(np.pi/2 - np.arcsin(na1/n))
    theta1 = np.pi/2 - np.arcsin(na1/n)
    theta2 = -(np.pi/2 - np.arcsin(na2/n))
    total_photons = 500
    dose_ratio = 10**param[0]
    dose2 = total_photons/(1 + dose_ratio)
    dose1 = total_photons - dose2
    
    exp = multiframe.MultiFrameMicroscope(ill_thetas=[theta1, theta2], det_thetas=[theta2, theta1],
                                      ill_nas=[na1, na2], det_nas=[na2, na1],
                                      ill_types=2*['sheet'], det_types=2*['lens'],
                                      colors=['(1,0,0)', '(0,0,1)'], n_frames=4,
                                      n_pts=n_pts, max_photons=[dose1, dose2], n_samp=1.33)

    exp.calc_estimation_stats()
    data = exp.sa_uncert
    med = np.median(data)
    return med, np.median(np.abs(data - med))

med = []
mad = []
for i, pt in enumerate(pts.T):
    print('Calculating microscope '+str(i+1)+'/'+str(pts.shape[1]))
    x = calc_stats(pt)
    med.append(x[0])
    mad.append(x[1])

# Plot 2D regions
def plot_2d_regions(ax, cax, pts, data, special_pt=(-1,-1)):
    
    # Set y ticks
    from matplotlib.ticker import FuncFormatter, FixedLocator
    ax.set_yticks([0, 30, 60, 90])
    ax.set_yticklabels(['0${}^{\circ}$', '30${}^{\circ}$', '60${}^{\circ}$', '90${}^{\circ}$'])
    ax2 = ax.twinx()
    ax2.set_yticks([0, 30, 60, 90])
    ax2.set_yticklabels([' 90${}^{\circ}$', '60${}^{\circ}$', '30${}^{\circ}$', '0${}^{\circ}$'])

    ax2.tick_params(axis='both', labelsize=14)

    xticks = [-2, -1, 0, 1, 2]
    ax.set_xticks(xticks)
    ax.set_xticklabels(['$10^{'+str(tick)+'}$' for tick in xticks])


    # Annotation
    def my_annotate(ax, annotation, xy, fontsize=9, rotation=0):
        ax.annotate(annotation, xy=(0,0), xytext=xy, textcoords='axes fraction',
                    va='center', ha='center', fontsize=fontsize,
                    annotation_clip=False, rotation=rotation, zorder=13)

    my_annotate(ax, 'Exposure${}_1$/Exposure${}_2$', (0.5, -0.13), fontsize=14)
    my_annotate(ax, '$\\alpha_{1}$, Objective 1 Half Angle', (-0.18, 0.5), fontsize=14, rotation=90)
    my_annotate(ax, '$\\alpha_{2}$, Objective 2 Half Angle', (1.17, 0.5), fontsize=14, rotation=-90)

    # Calculate colors
    color_map='coolwarm'
    color_norm='log'
    color_min=1e-4
    color_max=1e1
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

    ax.plot(special_pt[0], special_pt[1], 'kx', markersize=5, zorder=15)

    # Plot patches
    width = dose_width/(n_grid_pts)
    height = alpha_max/(n_grid_pts)
    for i, (pt, c) in enumerate(zip(pts_list, colors)):
        ax.add_patch(patches.Rectangle((pt[0] - width/2, pt[1] - height/2), width, height, facecolor=c, edgecolor=c))
    fig.colorbar(sc, cax=cax, orientation='vertical')

    # Mask around lines
    ax.set(xlim=[-2, 2], ylim=[0, 90])

# Plot first two columns
alpha1 = 60
na1 = 1.33*np.sin(np.deg2rad(alpha1))
na2 = n*np.sin(np.pi/2 - np.arcsin(na1/n))
theta1 = np.pi/2 - np.arcsin(na1/n)
theta2 = -(np.pi/2 - np.arcsin(na2/n))
dose_rat_plot = 0.6985
dose_ratio = 10**dose_rat_plot
total_photons = 500
dose2 = total_photons/(1 + dose_ratio)
dose1 = total_photons - dose2

exp = multiframe.MultiFrameMicroscope(ill_thetas=[theta1, theta2], det_thetas=[theta2, theta1],
                                      ill_nas=[na1, na2], det_nas=[na2, na1],
                                      ill_types=2*['sheet'], det_types=2*['lens'],
                                      colors=['(1,0,0)', '(0,0,1)'], n_frames=4,
                                      n_pts=n_pts_sphere, max_photons=[dose1, dose2], n_samp=1.33)
exp.calc_estimation_stats()

# Make scene string
# alpha1 label
scene_string = exp.scene_string()
line_string = "draw(O--expi(theta, 0));\n"
line_string = line_string.replace('theta', str(np.deg2rad(90-alpha1)))
scene_string += line_string
line_string = "draw(O--X);\n"
scene_string += line_string
arc_string = 'draw(L=Label("$\\alpha_1$", align=1.5*W), arc(O, 0.1*expi(theta, 0), 0.1*X, normal=Y));\n'
arc_string = arc_string.replace('theta', str(np.deg2rad(90-alpha1)))
scene_string += arc_string

# alpha2 label
line_string = "draw(O--expi(theta, 0));\n"
line_string = line_string.replace('theta', str(np.deg2rad(-alpha1)))
scene_string += line_string
line_string = "draw(O--(-X));\n"
scene_string += line_string
arc_string = 'draw(L=Label("$\\alpha_2$", align=7*NE), arc(O, -0.1*X, 0.1*expi(theta, 0), normal=Y));\n'
arc_string = arc_string.replace('theta', str(np.deg2rad(-alpha1)))
scene_string += arc_string

# Perpendicular 
line_string = "draw(0.075*expi(theta2, 0)--(0.075*expi(theta2, 0)+0.075*expi(theta1, 0)));\n"
line_string = line_string.replace('theta1', str(np.deg2rad(90-alpha1)))
line_string = line_string.replace('theta2', str(np.deg2rad(-alpha1)))
scene_string += line_string
line_string = "draw(0.075*expi(theta1, 0)--(0.075*expi(theta1, 0)+0.075*expi(theta2, 0)));\n"
line_string = line_string.replace('theta1', str(np.deg2rad(90-alpha1)))
line_string = line_string.replace('theta2', str(np.deg2rad(-alpha1)))
scene_string += line_string

util.draw_scene(scene_string, my_ax=ax0, dpi=dpi)
util.draw_scene(exp.ellipse_string(n_pts=250), my_ax=ax1, dpi=dpi)
util.plot_sphere(directions=exp.directions, data=exp.sa_uncert,
                 color_norm='log', linthresh=1e-4,
                 color_max=1e1, color_min=1e-4,
                 my_ax=ax2, my_cax=cax2)
    
# Plots last two columns
plot_2d_regions(ax3, cax3, pts, med, special_pt=(dose_rat_plot, alpha1))
plot_2d_regions(ax4, cax4, pts, mad, special_pt=(dose_rat_plot, alpha1))
    
# Label axes and save
print('Saving final figure.')    
fig.savefig('../paper/asymmetric-double.pdf', dpi=250)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
