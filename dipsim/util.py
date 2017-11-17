from dipsim import visuals
import numpy as np
import subprocess
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import make_axes_locatable
import vispy
from skimage import io

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
    
def tp2xyz(theta, phi):
    # Convert spherical to cartesian
    return np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)

def xyz2tp(x, y, z):
    # Convert cartesian to spherical
    arccos_arg = z/np.sqrt(x**2 + y**2 + z**2)
    if np.isclose(arccos_arg, 1.0): # Avoid arccos floating point issues
        arccos_arg = 1.0
    elif np.isclose(arccos_arg, -1.0):
        arccos_arg = -1.0
    return np.arccos(arccos_arg), np.arctan2(y, x)

def axis_angle(theta1, phi1, theta2, phi2):
    # Find angle between axes in spherical coordinates
    arccos_arg = np.abs(np.dot(tp2xyz(theta1, phi1), tp2xyz(theta2, phi2)))
    if np.isclose(arccos_arg, 1.0, atol=1e-5): # Avoid arccos floating point issues
        arccos_arg = 1.0
    elif np.isclose(arccos_arg, -1.0):
        arccos_arg = -1.0

    return np.arccos(arccos_arg)

# Plotting functions
def fibonacci_sphere(n):
    # Returns "equally" spaced points on a unit sphere in spherical coordinates.
    # http://stackoverflow.com/a/26127012/5854689
    z = np.linspace(1 - 1/n, -1 + 1/n, num=n) 
    theta = np.arccos(z)
    phi = np.mod((np.pi*(3.0 - np.sqrt(5.0)))*np.arange(n), 2*np.pi) - np.pi
    return np.vstack((theta, phi)).T

def sphere_profile(n):
    # Returns equally spaced points on a profile around the unit sphere in
    # spherical coordinates. Only on the equator for now. 
    z = np.linspace(1 - 1/n, -1 + 1/n, num=n)
    theta = n*[np.pi/2]
    phi = np.linspace(-np.pi, np.pi, num=n)
    return np.vstack((theta, phi)).T


def plot_sphere(filename=None, directions=None, data=None, interact=False,
                color_norm='linear', color_min=0, color_max=None,
                gamma=0.25, color_map='coolwarm', linthresh=1e-3,
                my_ax=None, my_cax=None, dpi=500, vis_px=1000):

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
        if filename is not None:
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

def tiff2array(filename, x=0, y=0, z=0, width=None, height=None, slices=None):
    # Opens .tiff and returns array starting at x, y, z with width, height, and
    # slices dimensions. "None" means return the whole dimension.
    im = io.imread(filename)
    shape = im.shape
    x_min = x
    if width is None:
        x_max = shape[2]
    else:
        x_max = x + width
    y_min = y
    if height is None:
        y_max = shape[1]
    else:
        y_max = y + height
    z_min = z
    if slices is None:
        z_max = shape[0]
    else:
        z_max = z + slices
    im = im[z_min:z_max, y_min:y_max, x_min:x_max]
    if slices == 1:
        return im.squeeze()
    else:
        return im

def plot_array(input_files, output_file, row_labels=[], col_labels=[],
               line_start=(0,0), line_end=(0,0),
               roi_upper_left=(0, 0), roi_wh=(0,0), plot_roi=False, zslice=0):
    # Create axes
    shape = input_files.shape
    inches = 3
    fig, axs = plt.subplots(2*shape[0], shape[1], figsize=(3.5*inches, 3*inches), gridspec_kw={'height_ratios':[1, 0.2, 1, 0.2], 'hspace':0.1, 'wspace':0.1})
    main_axs = axs[::2,:]
    line_axs = axs[1::2,:]

    # Load data and plot
    for (input_file, main_ax, line_ax) in zip(input_files.flatten(), main_axs.flatten(), line_axs.flatten()):
        # Load data
        if plot_roi:
            im = tiff2array(input_file, x=roi_upper_left[0], y=roi_upper_left[1], z=zslice, 
                            width=roi_wh[0], height=roi_wh[1], slices=1)
        else:
            im = tiff2array(input_file, z=zslice, slices=1)
            
        # Plot image
        main_ax.imshow(im, interpolation=None)
        main_ax.set_axis_off()
        line_ax.set_axis_off()

        # Plot roi box on image
        if not plot_roi:
            rx0, ry0 = roi_upper_left
            rx1, ry1 = rx0 + roi_wh[0], ry0 + roi_wh[1],
            main_ax.plot([rx0, rx1, rx1, rx0, rx0], [ry0, ry0, ry1, ry1, ry0], 'g-')
        
        # Plot line on image
        x0, y0 = line_start
        x1, y1 = line_end
        main_ax.plot([x0, x1], [y0, y1], 'r-')
        main_ax.plot([x0], [y0], 'r.', markersize=6)

        # Extract line profile
        num = 1000
        x, y = np.linspace(x0, x1, num), np.linspace(y0, y1, num)
        zi = []
        for (xi, yi) in zip(x, y):
            zi.append(im[int(yi), int(xi)])

        # Plot profile
        line_ax.plot(zi, '-r')
        line_ax.plot(zi[0], '.r', markersize=6)
        line_ax.set_ylim(0, 2000)
        line_ax.set_xlim(-100, 1100)
        print(input_file, main_ax)

    # Label rows and columns
    for i, label in enumerate(col_labels):
        axs[0, i].annotate(label, xy=(0,0), xytext=(0.5, 1.1), textcoords='axes fraction', va='center', ha='center', fontsize=12, annotation_clip=False)
    for i, label in enumerate(row_labels):
        axs[2*i, 0].annotate(label, xy=(0,0), xytext=(-0.1, 0.5), textcoords='axes fraction', va='center', ha='center', fontsize=12, annotation_clip=False, rotation=90)

    fig.savefig(output_file, dpi=200)
