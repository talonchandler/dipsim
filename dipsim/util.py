from dipsim import visuals
import numpy as np
import subprocess
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import make_axes_locatable
import vispy

# Coordinate conversion functions. 
def theta_prime(theta, phi, psi):
    return np.arccos(np.sin(psi)*np.cos(phi)*np.sin(theta) + np.cos(psi)*np.cos(theta))

def phi_prime(theta, phi, psi):
    num = np.cos(psi)*np.cos(phi)*np.sin(theta) - np.sin(psi)*np.cos(theta)
    den = np.sqrt(1 - (np.sin(psi)*np.cos(phi)*np.sin(theta) + np.cos(psi)*np.cos(theta))**2)
    if phi < np.pi and phi > 0:
        return np.arccos(num/den)
    elif phi < 0 and phi > -np.pi:
        return -np.arccos(num/den)
    else:
        return 0
    
def tp2xyz(tp, R=np.eye(3,3)):
    xyz = [np.sin(tp[0])*np.cos(tp[1]), np.sin(tp[0])*np.sin(tp[1]), np.cos(tp[0])]
    return np.dot(R, xyz)

def xyz2tp(xyz, R=np.eye(3,3)):
    xyz_prime = np.dot(np.linalg.inv(R), xyz)
    theta = np.arccos(xyz_prime[2]/np.linalg.norm(xyz_prime)),
    phi = np.arctan2(xyz_prime[1], xyz_prime[0])
    return np.array([theta, phi])

# Plotting functions
def fibonacci_sphere(n):
    # Returns "equally" spaced points on a unit sphere in spherical coordinates.
    # http://stackoverflow.com/a/26127012/5854689
    z = np.linspace(1 - 1/n, -1 + 1/n, num=n) 
    theta = np.arccos(z)
    phi = np.mod((np.pi*(3.0 - np.sqrt(5.0)))*np.arange(n), 2*np.pi) - np.pi
    return np.vstack((theta, phi)).T

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

def coeff_of_variation(data):
    return np.std(data)/np.mean(data)

def plot_histogram(data, ax):
    bins = np.array([10**x for x in np.arange(-4, 2, 0.1)])
    hist, bin_edges = np.histogram(data, bins=bins)
    ax.step(bin_edges[:-1], hist/len(data), '-k')
    ax.set_xscale('log')
    ax.set_yscale('log')    
    ax.set_xlim([np.min(data/10), 4*np.pi])
    ax.set_ylim([1e-3, 1])
    ax.set_xlabel(r'$\sigma_{\Omega}$', fontsize=18)
    cv = coeff_of_variation(data)
    cv_str = '$c_v =' + '{:.1e}'.format(cv).replace('e', '\\times 10^{') + '}$'
    ax.annotate(cv_str, xy=(0,0), xytext=(0.95, 0.95), textcoords='axes fraction',
                           va='center', ha='right', fontsize=16, annotation_clip=False)
    return cv
            
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

def draw_scene(scene_string, filename='out.png', my_ax=None, dpi=500,
               save_file=False, chop=True):
                
    asy_string = """
    import three;
    import graph3;
    settings.outformat = "pdf";
    settings.prc = true;
    settings.embed= true;
    settings.render=16;

    size(6cm,6cm);
    currentprojection = orthographic(1, 1, 1);

    void circle(real Theta, real Alpha, bool dash, triple color) {
      triple normal = expi(Theta, 0);
      real h = 1 - sqrt(2 - 2*cos(Alpha) - sin(Alpha)^2);
      real radius = sin(Alpha);
      path3 mycircle = circle(c=h*normal, r=radius, normal=normal);
      if (dash) {
        draw(mycircle, p=linetype(new real[] {8,8}, offset=xpart(color))+rgb(xpart(color), ypart(color), zpart(color)));
      } else {
        draw(mycircle, p=rgb(xpart(color), ypart(color), zpart(color)));
      }
    }

    void ellipse(real Theta, real Phi, real a, real b, real theta, bool dash, triple color) {
      triple normal = expi(Theta, Phi);
      real a_scaled = a/max(a, b);
      real b_scaled = b/max(a, b);      
      path3 mycircle = rotate(degrees(Phi), Z)*rotate(degrees(Theta), Y)*shift(Z)*rotate(degrees(theta), Z)*scale(a_scaled, b_scaled, 1)*circle(c=O, r=0.05, normal=Z);
      if (dash) {
        draw(mycircle, p=linetype(new real[] {8,8}, offset=xpart(color))+rgb(xpart(color), ypart(color), zpart(color)));
      } else {
        draw(mycircle, p=rgb(xpart(color), ypart(color), zpart(color)));
      }
    }

    void mydot(real Theta, triple color) {
      triple normal = expi(Theta, 0);
      dot(normal, p=rgb(xpart(color), ypart(color), zpart(color)));
    }

    void arrow(real Theta, real Phi_Pol, triple color, bool dash) {
      if (dash) {
        draw(rotate(Theta, Y)*rotate(Phi_Pol, Z)*(Z--(Z+0.2*X)), p=linetype(new real[] {4,4}, offset=xpart(color))+rgb(xpart(color), ypart(color), zpart(color)), arrow=Arrow3(emissive(rgb(xpart(color), ypart(color), zpart(color)))));
        draw(rotate(Theta, Y)*rotate(Phi_Pol, Z)*(Z--(Z-0.2*X)), p=linetype(new real[] {4,4}, offset=xpart(color))+rgb(xpart(color), ypart(color), zpart(color)), arrow=Arrow3(emissive(rgb(xpart(color), ypart(color), zpart(color)))));
      } else {
        draw(rotate(Theta, Y)*rotate(Phi_Pol, Z)*(Z--(Z+0.2*X)), p=rgb(xpart(color), ypart(color), zpart(color)), arrow=Arrow3(emissive(rgb(xpart(color), ypart(color), zpart(color)))));
        draw(rotate(Theta, Y)*rotate(Phi_Pol, Z)*(Z--(Z-0.2*X)), p=rgb(xpart(color), ypart(color), zpart(color)), arrow=Arrow3(emissive(rgb(xpart(color), ypart(color), zpart(color)))));
      }
    }

    void watson(real Theta, real Phi, real kappa, real x, real y, real z) {
     int n_phi = 10;
     int n_theta = 10;

     real max_radius = 0;
     if(kappa > 0){
       max_radius = exp(kappa);
     }
     else{
       max_radius = 1.0;
     }

     for(int i=0; i <= n_theta; ++i) {
       real Theta_i = pi*i/n_theta;
       real weight = exp(kappa*(cos(Theta_i)**2))/max_radius;     
       path3 mycircle = circle(c=Z*weight*cos(Theta_i), r=weight*sin(Theta_i));
       draw(shift((x, y, z))*rotate(angle=degrees(Phi), u=O, v=Z)*rotate(angle=degrees(Theta), u=O, v=Y)*mycircle);
     }

     triple f(real t) {
       real weight = exp(kappa*(cos(t)**2))/max_radius;
       return (0, weight*sin(t), weight*cos(t));
     }
     path3 phi_path = graph(f, 0, 2pi, operator ..);

     for(int i=0; i <= n_phi; ++i) {
       real Phi_i = 2*pi*i/n_theta;
       draw(shift((x, y, z))*rotate(angle=degrees(Phi), u=O, v=Z)*rotate(angle=degrees(Theta), u=O, v=Y)*rotate(angle=degrees(Phi_i), u=(0,0,0), v=(0,0,1))*phi_path);
     }
    }
    real len = 10;
    draw((-len,-len)--(len,-len)--(len,len)--(-len,len)--(-len,-len), white);
    """

    asy_string += scene_string
    asy_string += "dot(Z);shipout(scale(4.0)*currentpicture.fit());"

    text_file = open("temp.asy", "w")
    text_file.write(asy_string)
    text_file.close()

    subprocess.call(['asy', 'temp.asy'])
    subprocess.call(['convert', '-density', str(dpi), '-units', 'PixelsPerInch', 'temp.pdf', 'temp.png'])

    im = mpimg.imread('temp.png')

    # Chop top of im to make it square and fix asy error
    if chop:
        im = im[int(im.shape[1]*0.075):,:,:]
    
    f = plt.figure(figsize=(5, 5), frameon=False)
    local_ax = plt.axes([0, 0, 1, 1]) # x, y, width, height
    if my_ax == None:
        my_ax = local_ax

    for ax in [local_ax, my_ax]:
        draw_axis(ax)
        ax.spines['right'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])

        # Plot
        ax.imshow(im, interpolation='none')

    # Save
    if save_file:
        f.savefig(filename, dpi=dpi)

    subprocess.call(['rm', 'temp.asy', 'temp.pdf', 'temp.png'])
    return ax
