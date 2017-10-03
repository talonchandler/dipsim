from dipsim import visuals
import numpy as np
import subprocess
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import make_axes_locatable
import vispy

def plot_sphere(filename=None, directions=None, data=None, interact=False,
                color_norm='linear', color_min=0, color_max=None,
                gamma=0.25, color_map='coolwarm', linthresh=1e-3,
                my_ax=None, my_cax=None, dpi=500, vis_px=1000,
                save_file=False):

    # Setup viewing window
    vispy.use('glfw')
    canvas = vispy.scene.SceneCanvas(keys='interactive', bgcolor='white',
                               size=(vis_px, vis_px), show=interact)
    my_cam = vispy.scene.cameras.turntable.TurntableCamera(fov=0, azimuth=135,
                                                           scale_factor=2.05)
    view = canvas.central_widget.add_view(camera=my_cam)

    # Plot dots
    dots = vispy.scene.visuals.Markers(parent=view.scene)
    dots.antialias = False
    dots.set_data(pos=np.array([[1.01,0,0],[0,1.01,0],[0,0,1.01]]),
                  edge_color='black', face_color='black', size=vis_px/50)

    if color_max == None:
        color_max = np.max(data)
    
    # Calculate colors
    if color_norm == 'linear':
        norm = matplotlib.colors.Normalize(vmin=color_min, vmax=color_max)
    elif color_norm == 'log':
        norm = matplotlib.colors.LogNorm(vmin=color_min, vmax=color_max)
    elif color_norm == 'linlog':
        norm = matplotlib.colors.SymLogNorm(linthresh=linthresh, vmin=-color_max, vmax=color_max)
    elif color_norm == 'power':
        norm = matplotlib.colors.PowerNorm(gamma=gamma, vmin=data.min(), vmax=data.max())

    norm_data = norm(data).data
    norm_data = np.expand_dims(norm_data, 1)    
    cmap = matplotlib.cm.get_cmap(color_map)
    colors = np.apply_along_axis(cmap, 1, norm_data)
    
    # Plot sphere
    sphere = visuals.MySphere(parent=view.scene, radius=1.0,
                              directions=directions, colors=colors)
    
    # Display or save
    if interact:
        visuals.MyXYZAxis(parent=view.scene, origin=[0,1.3,-0.3], length=0.2)
        canvas.app.run()
    else:
        im = canvas.render()
        f = plt.figure(figsize=(5,5))
        local_ax = plt.axes([0.0, 0, 0.9, 0.9]) # x, y, width, height
        local_cax = plt.axes([0.95, 0.01, 0.025, 0.86])

        if my_ax == None:
            my_ax = local_ax
        if my_cax == None:
            my_cax = local_cax
        
        for (ax, cax) in [(local_ax, local_cax), (my_ax, my_cax)]:
            ax.axis('off')
            draw_axis(ax)
            cmap = ax.imshow(im, interpolation='none', cmap=color_map, norm=norm)
            f.colorbar(cmap, cax=cax, orientation='vertical')

        # Save
        if save_file:
            f.savefig(filename, dpi=dpi)
            
def draw_axis(ax, x=0.925, y=0.1):
    length=0.1
    center=np.array((x, y))
    angles = np.pi/2 + np.array((0, 2*np.pi/3, 4*np.pi/3))
    labels = ('$z$', '$x$', '$y$')
    for angle, label in zip(angles, labels):
        end = center + np.array((length*np.cos(angle), length*np.sin(angle)))
        ax.annotate(label,ha='center', va='center', xy=center, xytext=end,
                    xycoords='axes fraction', 
                    arrowprops=dict(arrowstyle="<-", shrinkA=0, shrinkB=0),
                    fontsize=14)

def generate_caxs(axs):
    caxs = []
    for ax in axs.flatten():
        divider = make_axes_locatable(ax)
        caxs.append(divider.append_axes("right", size="5%", pad=0.15))
    return np.array(caxs).reshape(axs.shape)
        
