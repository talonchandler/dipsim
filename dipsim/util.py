from dipsim import visuals
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import vispy

def normalize(x):
    """ 
    Returns a normalized vector. Returns zero vector if input is zero.
    """
    len_x = np.linalg.norm(x)
    if len_x == 0:
        return x
    else:
        return x/len_x

def rot_mat(theta, u):
    """
    Returns the rotation matrix that performs a right handed rotation by 
    angle theta about the vector u.

    Reference: https://en.wikipedia.org/wiki/Rodrigues'_rotation_formula
    """
    u = u/np.linalg.norm(u)
    K = np.array([[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]])
    return np.identity(3) + K*np.sin(theta) + np.dot(K, K)*(1 - np.cos(theta))

def rot_map(u, v=np.array([0,0,1])):
    """
    Returns the rotation matrix that aligns v with u.
    """
    if np.array_equal(u, v):
        return np.eye(3)
    else:
        return rot_mat(np.arccos(np.dot(v, u)), np.cross(v, u))

def orthonormal_basis(v0):
    """
    Returns two orthonormal vectors that are orthogonal to v0.
    """
    if np.dot(v0, [0, 0, 1]) == 0:
        v1 = np.array([0, 0, 1])                
    else:
        v1 = np.array([1, 0, -v0[0]/v0[2]])
    v1 = v1/np.linalg.norm(v1)
    v2 = np.cross(v1, v0)
    v2 = v2/np.linalg.norm(v2)
    return v1, v2

def my_orthonormal_basis(v0):
    """
    Returns two orthonormal vectors that are orthogonal to v0.
    """
    if np.dot(v0, [0, 0, 1]) == 1:
        v1 = np.array([1, 0, 0])                
    else:
        v1 = np.array([-v0[1], v0[0], 0])
    v1 = v1/np.linalg.norm(v1)
    v2 = np.cross(v0, v1)
    v2 = v2/np.linalg.norm(v2)
    return v1, v2

def fibonacci_sphere(n):
    # Returns "equally" spaced points on a unit sphere in spherical coordinates.
    # http://stackoverflow.com/a/26127012/5854689
    z = np.linspace(1 - 1/n, -1 + 1/n, num=n) 
    theta = np.arccos(z)
    phi = np.mod((np.pi*(3.0 - np.sqrt(5.0)))*np.arange(n), 2*np.pi) - np.pi
    return np.vstack((theta, phi)).T

# Three coordinate conversion functions. Use R to use theta-phi coordinates that
# are measured from axes other than the typical z and x axes.
def tp2xyz(tp, R=np.eye(3,3)):
    xyz = [np.sin(tp[0])*np.cos(tp[1]), np.sin(tp[0])*np.sin(tp[1]), np.cos(tp[0])]
    return np.dot(R, xyz)

def xyz2tp(xyz, R=np.eye(3,3)):
    xyz_prime = np.dot(np.linalg.inv(R), xyz)
    theta = np.arccos(xyz_prime[2]/np.linalg.norm(xyz_prime)),
    phi = np.arctan2(xyz_prime[1], xyz_prime[0])
    # if phi < 0:
    #     phi += 2*np.pi
    return np.array([theta, phi])

def tp2tp_prime(tp, R=np.eye(3,3)):
    if np.array_equal(R, np.eye(3,3)):
        return tp
    else:
        xyz = tp2xyz(tp)
        return xyz2tp(xyz, R)

# Plotting functions
def plot_sphere(filename=None, directions=None, data=None, interact=False,
                color_norm='linear', color_min=0, color_max=4*np.pi,
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

def dispersion_index(data):
    return np.var(data)/np.mean(data)
            
def plot_histogram(data, ax):
    bins = np.array([10**x for x in np.arange(-4, 2, 0.1)])
    hist, bin_edges = np.histogram(data, bins=bins)
    ax.step(bin_edges[:-1], hist/len(data), '-k')
    ax.set_xscale('log')
    ax.set_yscale('log')    
    ax.set_xlim([np.min(data/10), 4*np.pi])
    ax.set_ylim([1/(2*len(data)), 1])
    ax.set_xlabel(r'$\sigma_{\Omega}$', fontsize=18)
    d = dispersion_index(data)
    d_str = '$D =' + '{:.1e}'.format(d).replace('e', '\\times 10^{') + '}$'
    ax.annotate(d_str, xy=(0,0), xytext=(0.95, 0.95), textcoords='axes fraction',
                           va='center', ha='right', fontsize=16, annotation_clip=False)
    return d
            
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

def label_rows_and_cols(axs, row_labels='', col_labels='',
                        row_pos=(-0.1, 0.5), col_pos=(0.5, 1.1)):
    for i, label in enumerate(row_labels):
        axs[i][0].annotate(label, xy=(0,0), xytext=row_pos, textcoords='axes fraction',
                           va='center', ha='center', rotation=90, fontsize=18, annotation_clip=False)
    for i, label in enumerate(col_labels):
        axs[0][i].annotate(label, xy=(0,0), xytext=col_pos, textcoords='axes fraction',
                           va='center', ha='center', fontsize=18, annotation_clip=False)

def generate_caxs(axs):
    caxs = []
    for ax in axs.flatten():
        divider = make_axes_locatable(ax)
        caxs.append(divider.append_axes("right", size="5%", pad=0.15))
    return np.array(caxs).reshape(axs.shape)

def a3to6(vec3):
    # Converts a gjv from three components to six components.
    return np.array([np.conj(vec3[0])*vec3[0],
                     np.conj(vec3[1])*vec3[1],
                     np.conj(vec3[2])*vec3[2],
                     2*np.real(np.conj(vec3[0])*vec3[1]),
                     2*np.real(np.conj(vec3[0])*vec3[2]),
                     2*np.real(np.conj(vec3[1])*vec3[2])])

def mu3to6(vec3):
    # Converts a dipole vector from three components to six components.
    return np.array([np.conj(vec3[0])*vec3[0],
                     np.conj(vec3[1])*vec3[1],
                     np.conj(vec3[2])*vec3[2],
                     np.real(np.conj(vec3[0])*vec3[1]),
                     np.real(np.conj(vec3[0])*vec3[2]),
                     np.real(np.conj(vec3[1])*vec3[2])])
