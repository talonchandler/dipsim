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

    void dip_arrow(real X, real Y, real Theta, real Phi, real length, triple color) {
      //draw(-expi(Theta, Phi)--expi(Theta, Phi));
      draw(shift((X, Y, 0))*(O--length*expi(Theta, Phi)), p=rgb(xpart(color), ypart(color), zpart(color)), arrow=Arrow3(emissive(rgb(xpart(color), ypart(color), zpart(color)))));
      draw(shift((X, Y, 0))*(O--(-length*expi(Theta, Phi))), p=rgb(xpart(color), ypart(color), zpart(color)), arrow=Arrow3(emissive(rgb(xpart(color), ypart(color), zpart(color)))));
      dot((X,Y,0), p=rgb(xpart(color), ypart(color), zpart(color)));
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
    draw(scale(128, 256, 0)*unitplane, surfacepen=white+opacity(0.5));
    defaultpen(fontsize(8pt));
    draw((0, -10, 0)--(50, -10, 0), L=Label("$50$", position=MidPoint, align=NW));
    dotfactor=2;
    """

    asy_string += scene_string
    asy_string += "dot(O);shipout(scale(4.0)*currentpicture.fit());"

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
        print("Saving: "+filename)
        f.savefig(filename, dpi=dpi)

    subprocess.call(['rm', 'temp.asy', 'temp.pdf', 'temp.png'])
    return ax
