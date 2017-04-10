import numpy as np
import matplotlib.pyplot as plt
import functools
import vispy
from dipsim import util, fluorophore, visuals
from vispy.visuals.transforms import (STTransform, LogTransform,
                                      MatrixTransform, PolarTransform)

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

    def calc_intensity(self, fluorophores):
        ill = self.illuminator
        det = self.detector
        I = 0
        for f in fluorophores:
            excite = self.illuminator.calc_excitation_efficiency(flu)
            collect = self.detector.calc_collection_efficiency(flu)
            I += excite*collect
        return I*self.max_photons

    def calc_sensitivity(self, args):
        theta = args[0]
        phi = args[1]

        flu_dir = np.array([np.sin(theta)*np.cos(phi),
                           np.sin(theta)*np.sin(phi),
                           np.cos(theta)])

        flu = fluorophore.Fluorophore(position=np.array([0, 0, 0]),
                                      mu_abs=flu_dir,
                                      mu_em=flu_dir)

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
        theta = args[0]
        phi = args[1]

        flu_dir = np.array([np.sin(theta)*np.cos(phi),
                           np.sin(theta)*np.sin(phi),
                           np.cos(theta)])

        flu = fluorophore.Fluorophore(position=np.array([0, 0, 0]),
                                      mu_abs=flu_dir,
                                      mu_em=flu_dir)

        I = self.illuminator.calc_excitation_efficiency(flu)
        return I
    
    def plot_excitation_efficiency(self, filename='out.png',
                                                 n=50, **kwargs):
        directions = util.fibonacci_sphere(n)
        print('Generating data for microscope: '+filename)
        I = np.apply_along_axis(self.calc_excitation_efficiency,
                                      1, directions)
        print('Plotting data for microscope: '+filename)
        util.plot_sphere(filename, directions=directions, data=I, **kwargs)

    def calc_collection_efficiency(self, args):
        theta = args[0]
        phi = args[1]

        flu_dir = np.array([np.sin(theta)*np.cos(phi),
                           np.sin(theta)*np.sin(phi),
                           np.cos(theta)])

        flu = fluorophore.Fluorophore(position=np.array([0, 0, 0]),
                                      mu_abs=flu_dir,
                                      mu_em=flu_dir)

        I = self.detector.calc_collection_efficiency(flu)
        return I
    
    def plot_collection_efficiency(self, filename='out.png',
                                                 n=50, **kwargs):
        directions = util.fibonacci_sphere(n)
        print('Generating data for microscope: '+filename)
        I = np.apply_along_axis(self.calc_collection_efficiency,
                                      1, directions)
        print('Plotting data for microscope: '+filename)
        util.plot_sphere(filename, directions=directions, data=I, **kwargs)
        
    def draw_scene(self, filename='out.png', interact=False, my_ax=None, dpi=500,
                   vis_px=2000, save_file=False, pol_dirs=None):
        vispy.use('glfw')
        vis = vispy.scene.visuals
        
        # Setup viewing window
        canvas = vispy.scene.SceneCanvas(keys='interactive', bgcolor='white',
                                         size=(vis_px, vis_px), show=interact)
        cam = vispy.scene.cameras.turntable.TurntableCamera
        my_cam = cam(fov=0, azimuth=135, scale_factor=23,
                     center=(0, 0, self.illuminator.f+3))
        view = canvas.central_widget.add_view(camera=my_cam)

        # Plot dipole
        dip = visuals.MyArrow(parent=view.scene, length=2)
        m = MatrixTransform()
        m.rotate(angle=-45, axis=(1,1,0))
        dip.transform = m
        
        # Plot illuminator
        i = self.illuminator
        
        # TODO: Generalize these transformations for oblique illuminations
        # TODO: Clean this up
        if pol_dirs == None:
            pol_dirs = [i.bfp_pol]
        for pol_dir in pol_dirs:
            for direction in [-1, 1]:
                pol = visuals.MyArrow(parent=view.scene, length=3)
                m = MatrixTransform()
                m.rotate(angle=direction*90, axis=(0,1,0))
                phi_angle = np.degrees(np.arctan2(pol_dir[1], pol_dir[0]))
                m.rotate(angle=phi_angle, axis=(0,0,1))
                m.translate((0,0,2*i.f+0.1))
                pol.transform = m
        
        circ = vis.Ellipse(parent=view.scene, center=(0,0), radius=i.bfp_rad,
                           color=(1, 1, 0, 0.2 + 0.8/i.bfp_rad))
        line = vis.Ellipse(parent=view.scene, radius=3.0,
                           color=(1, 1, 0, 0.2 + 0.8/i.bfp_rad))
        
        circ.transform = STTransform(translate=(0, 0, 2*i.f))
        
        # Plot lens
        lens2 = visuals.MyLens(parent=view.scene)
        lens2.transform = STTransform(scale=(10, 10, 10),
                                      translate=(0, 0, i.f))
        lens_out = vis.Ellipse(parent=view.scene, center=(0,0), radius=0.5,
                               color=None)
        lens_out.transform = STTransform(scale=(10, 10, 10),
                                         translate=(0, 0, i.f))

        # Display or save
        if interact:
            canvas.app.run()
        else:
            # Setup figure
            im = canvas.render()
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
