import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import functools
import subprocess
import vispy
from dipsim import util, fluorophore, visuals
from vispy.visuals.transforms import (STTransform, LogTransform,
                                      MatrixTransform, PolarTransform)

from vispy import gloo

class Microscope:
    """
    A Microscope represents an experiment that collects a single frame of 
    intensity data.  

    A Microscope is specified by its illumination path (an Illuminator object),
    and its detection path (a Detector object).
    """
    def __init__(self, illuminator, detector, max_photons):
        self.illuminator = illuminator
        self.detector = detector
        self.max_photons = max_photons

    def calc_intensity(self, args):
        return self.max_photons*self.calc_sensitivity(args)
    
    def calc_sensitivity(self, args):
        flu = fluorophore.Fluorophore(mu_abs=args, mu_em=args)
        excite = self.illuminator.calc_excitation_efficiency(flu)
        collect = self.detector.calc_collection_efficiency(flu)        
        return excite*collect

    def plot_sensitivity(self, filename='out.png', n=50, **kwargs):
        directions = util.fibonacci_sphere(n)
        print('Generating data for microscope: '+filename)
        I = np.apply_along_axis(self.calc_sensitivity, 1, directions)
        print('Plotting data for microscope: '+filename)
        util.plot_sphere(filename, directions=directions, data=I, **kwargs)

    def calc_excitation_efficiency(self, args):
        flu = fluorophore.Fluorophore(mu_abs=args, mu_em=args)
        I = self.illuminator.calc_excitation_efficiency(flu)
        return I
    
    def plot_excitation_efficiency(self, filename='out.png', n=50, **kwargs):
        directions = util.fibonacci_sphere(n)
        print('Generating data for microscope: '+filename)
        I = np.apply_along_axis(self.calc_excitation_efficiency,
                                      1, directions)
        print('Plotting data for microscope: '+filename)
        util.plot_sphere(filename, directions=directions, data=I, **kwargs)

    def calc_collection_efficiency(self, args):
        flu = fluorophore.Fluorophore(mu_abs=args, mu_em=args)        
        I = self.detector.calc_collection_efficiency(flu)
        return I
    
    def plot_collection_efficiency(self, filename='out.png', n=50, **kwargs):
        directions = util.fibonacci_sphere(n)
        print('Generating data for microscope: '+filename)
        I = np.apply_along_axis(self.calc_collection_efficiency,
                                      1, directions)
        print('Plotting data for microscope: '+filename)
        util.plot_sphere(filename, directions=directions, data=I, **kwargs)

    def draw_scene(self, filename='out.png', my_ax=None, dpi=500,
                    save_file=False, pol_dirs=None, dual_arm=False):

        asy_string = """
        import three;
        settings.outformat = "pdf";
        settings.prc = true;
        settings.embed= true;
        settings.render=16;

        size(6cm,0);
        currentprojection = orthographic(1, 1, 1);

        void circle(real Theta, real Alpha, bool dash, triple color) {
          triple normal = expi(Theta, 0);
          real h = 1 - sqrt(2 - 2*cos(Alpha) - sin(Alpha)^2);
          real radius = sin(Alpha);
          path3 mycircle = circle(c=h*normal, r=radius, normal=normal);
	  if (dash) {
	    draw(mycircle, p=dashed+rgb(xpart(color), ypart(color), zpart(color)));
	  } else {
	    draw(mycircle, p=rgb(xpart(color), ypart(color), zpart(color)));
	  }
        }

        void arrow(real Theta, real Phi_Pol, triple color) {
          draw(rotate(Theta, Y)*rotate(Phi_Pol, Z)*(Z--(Z+0.2*X)), p=rgb(xpart(color), ypart(color), zpart(color)), arrow=Arrow3(emissive(rgb(xpart(color), ypart(color), zpart(color)))));
          draw(rotate(Theta, Y)*rotate(Phi_Pol, Z)*(Z--(Z-0.2*X)), p=rgb(xpart(color), ypart(color), zpart(color)), arrow=Arrow3(emissive(rgb(xpart(color), ypart(color), zpart(color)))));
        }

        // Sphere
        draw(unitsphere, surfacepen=material(diffusepen=white+opacity(0.1), emissivepen=grey, specularpen=white));

        // Draw points on sphere
        dotfactor = 7;
        dot(X); 
        dot(Y); 
        dot(Z); 
        circle(0, pi/2, false, (0, 0, 0));
        """
        
        if dual_arm:
            dets = [self.illuminator, self.detector]
            ills = [self.detector, self.illuminator]
            colors = ['(1, 0, 0)', '(0, 0, 1)']
        else:
            dets = [self.detector]
            ills = [self.illuminator]
            colors = ['(1, 0, 0)']
            
        for idx, (det, ill, color) in enumerate(zip(dets, ills, colors)):

            # Plot illuminator
            illum_string = "circle(theta, alpha, false, color);\n"
            illum_string = illum_string.replace('theta', str(ill.theta_optical_axis))
            illum_string = illum_string.replace('alpha', str(ill.alpha))
            illum_string = illum_string.replace('color', str(color))        
            asy_string += illum_string

            # Plot detector
            detect_string = "circle(theta, alpha, true, color);\n"
            detect_string = detect_string.replace('theta', str(det.theta_optical_axis))
            detect_string = detect_string.replace('alpha', str(det.alpha+0.01))
            detect_string = detect_string.replace('color', color)        
            asy_string += detect_string

            # Plot polarizations
            if pol_dirs == None:
                pol_dirs = [self.illuminator.phi_pol]
            for pol_dir in pol_dirs:
                pol_string = "arrow(theta, phi_pol, color);\n"
                pol_string = pol_string.replace('theta', str(np.rad2deg(ill.theta_optical_axis)))
                pol_string = pol_string.replace('phi_pol', str(np.rad2deg(pol_dir)))
                pol_string = pol_string.replace('color', color)        
                asy_string += pol_string
        
        asy_string += "shipout(scale(4.0)*currentpicture.fit());"
        
        text_file = open("temp.asy", "w")
        text_file.write(asy_string)
        text_file.close()

        subprocess.call(['asy', 'temp.asy'])
        subprocess.call(['convert', '-density', str(dpi), '-units', 'PixelsPerInch', 'temp.pdf', 'temp.png'])
        #subprocess.call(['rm', 'temp.asy', 'temp.pdf', 'temp.png'])
        
        im = mpimg.imread('temp.png')
        f = plt.figure(figsize=(5, 5), frameon=False)
        local_ax = plt.axes([0, 0, 1, 1]) # x, y, width, height
        if my_ax == None:
            my_ax = local_ax

        for ax in [local_ax, my_ax]:
            util.draw_axis(ax)
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
        return ax
