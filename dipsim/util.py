from dipsim import visuals
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
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

def fibonacci_sphere(n):
    # Returns "equally" spaced points on a unit sphere in spherical coordinates. 
    pts = []
    offset = 2./n
    increment = np.pi * (3. - np.sqrt(5.));

    # TODO: Optimize this
    for i in range(n):
        y = ((i * offset) - 1) + (offset / 2);
        r = np.sqrt(1 - pow(y,2))

        phi = ((i) % n) * increment

        x = np.cos(phi) * r
        z = np.sin(phi) * r

        theta_out = np.arccos(z)        
        phi_out = np.arctan2(y, x)
        pts.append([theta_out, phi_out])

    return np.array(pts)

def tp2xyz(tp):
    return [np.sin(tp[0])*np.cos(tp[1]),
            np.sin(tp[0])*np.sin(tp[1]),
            np.cos(tp[0])]

def plot_sphere(filename, directions=None, data=None,
                interact=False, color_norm='linear', my_ax=None, my_cax=None,
                dpi=500, vis_px=1000):

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

    # Plot sphere
    sphere = visuals.MySphere(parent=view.scene, radius=1.0,
                              directions=directions, data=data,
                              color_norm=color_norm)
    
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
            
            # Colorbar
            if color_norm == 'linear':
                norm = matplotlib.colors.Normalize(vmin=data.min(), vmax=data.max())
                cmap = ax.imshow(im, interpolation='none', cmap='coolwarm', norm=norm)
                f.colorbar(cmap, cax=cax, orientation='vertical')
            elif color_norm == 'log':
                norm = matplotlib.colors.LogNorm(vmin=data.min(), vmax=data.max())
                cmap = ax.imshow(im, interpolation='none', cmap='coolwarm', norm=norm)
                cb = f.colorbar(cmap, cax=cax, orientation='vertical')
                cb.set_ticks([10**x for x in range(int(np.ceil(np.log10(data).min())), int(np.ceil(np.log10(data).max())))])
                
        # Save 
        f.savefig(filename, dpi=dpi)

def draw_axis(ax):
    length=0.1
    center=np.array((0.925, 0.1))
    angles = np.pi/2 + np.array((0, 2*np.pi/3, 4*np.pi/3))
    labels = ('$z$', '$x$', '$y$')
    for angle, label in zip(angles, labels):
        end = center + np.array((length*np.cos(angle), length*np.sin(angle)))
        ax.annotate(label,ha='center', va='center', xy=center, xytext=end,
                    xycoords='axes fraction', 
                    arrowprops=dict(arrowstyle="<-", shrinkA=0, shrinkB=0),
                    fontsize=14)
        
