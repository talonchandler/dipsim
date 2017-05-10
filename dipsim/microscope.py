import numpy as np
import matplotlib.pyplot as plt
import functools
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
        flu_dir = util.tp2xyz(args)
        flu = fluorophore.Fluorophore(mu_abs=flu_dir, mu_em=flu_dir)
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
        flu_dir = util.tp2xyz(args)
        flu = fluorophore.Fluorophore(mu_abs=flu_dir, mu_em=flu_dir)        
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
        flu_dir = util.tp2xyz(args)
        flu = fluorophore.Fluorophore(mu_abs=flu_dir, mu_em=flu_dir)        
        I = self.detector.calc_collection_efficiency(flu)
        return I
    
    def plot_collection_efficiency(self, filename='out.png', n=50, **kwargs):
        directions = util.fibonacci_sphere(n)
        print('Generating data for microscope: '+filename)
        I = np.apply_along_axis(self.calc_collection_efficiency,
                                      1, directions)
        print('Plotting data for microscope: '+filename)
        util.plot_sphere(filename, directions=directions, data=I, **kwargs)
        
    def draw_scene(self, filename='out.png', interact=False, my_ax=None, dpi=500,
                   vis_px=2000, save_file=False, pol_dirs=None, dual_arm=False):
        vispy.use('glfw')
        vis = vispy.scene.visuals
        i = self.illuminator
        k = np.array([0, 0, 1])
        
        # Setup viewing window
        canvas = vispy.scene.SceneCanvas(keys='interactive', bgcolor='white',
                                         size=(vis_px, vis_px), show=interact)
        cam = vispy.scene.cameras.turntable.TurntableCamera
        my_cam = cam(fov=0, azimuth=135, scale_factor=37, # large sf -> zoom out
                     center=(0, 0, self.illuminator.f))
        view = canvas.central_widget.add_view(camera=my_cam)

        if dual_arm:
            det_axes = [self.illuminator.optical_axis, self.detector.optical_axis]
            ill_axes = [self.detector.optical_axis, self.illuminator.optical_axis]
        else:
            det_axes = [self.detector.optical_axis]
            ill_axes = [self.illuminator.optical_axis]

        for idx, (det_axis, ill_axis) in enumerate(zip(det_axes, ill_axes)):
            # Plot illumination arm
            angle = -np.rad2deg(np.arccos(np.dot(ill_axis, k))),
            axis = np.cross(ill_axis, k)
            if np.linalg.norm(axis) == 0:
                axis = k
            
            # Plot illumination lens
            if idx == 0:
                lens = visuals.MyArrow(parent=view.scene, radius=i.bfp_rad, length=0.1,
                                       rows=1, cols=100, cone_length=0.01,
                                       color=(.8,.8,.8,1.0))
                m = MatrixTransform()
                m.translate((0, 0, i.f))
                m.rotate(angle=angle, axis=axis)
                lens.transform = m
            
        for idx, (det_axis, ill_axis) in enumerate(zip(det_axes, ill_axes)):
            # Detection arm
            angle = -np.rad2deg(np.arccos(np.dot(det_axis, k)))
            axis = np.cross(det_axis, k)            
            if np.linalg.norm(axis) == 0:
                axis = k

            # Axis
            ax = visuals.MyLine(parent=view.scene, length=i.f, radius=0.025)
            m = MatrixTransform()                        
            m.rotate(angle, axis)
            ax.transform = m

            # Collection lens
            if idx == 0:
                lens2 = vis.Ellipse(parent=view.scene, center=(0,0), radius=self.detector.det_rad,
                                   color=(.8,.8,.8,1.0))
                m = MatrixTransform()
                m.translate((0, 0, i.f))        
                m.rotate(angle, axis)
                lens2.transform = m

            # Plot detector
            if self.detector.det_type == '4pi':
                det = vis.Sphere(parent=view.scene, radius=3, color=(0.9,0.9,0.9,0.7))
            else:
                det = vis.Plane(parent=view.scene, width=2*self.detector.det_rad, height=2*self.detector.det_rad, color=(0.9,0.9,0.9,0.7))
                m = MatrixTransform()
                m.translate((0, 0, 2.0*i.f))          
                m.rotate(angle, axis)            
                det.transform = m

        for idx, (det_axis, ill_axis) in enumerate(zip(det_axes, ill_axes)):
            # Plot illumination arm
            angle = -np.rad2deg(np.arccos(np.dot(ill_axis, k))),
            axis = np.cross(ill_axis, k)
            if np.linalg.norm(axis) == 0:
                axis = k

            # Plot axis
            ax = visuals.MyLine(parent=view.scene, length=i.f, radius=0.025)
            m = MatrixTransform()
            m.rotate(angle=angle, axis=axis)
            ax.transform = m

            # # Plot illumination lens
            # if idx == 0:
            #     lens = visuals.MyArrow(parent=view.scene, radius=i.bfp_rad, length=0.1,
            #                            rows=1, cols=100, cone_length=0.01,
            #                            color=(.8,.8,.8,1.0))
            #     m = MatrixTransform()
            #     m.translate((0, 0, i.f))
            #     m.rotate(angle=angle, axis=axis)
            #     lens.transform = m

            # Plot illumination circle
            circ = visuals.MyArrow(parent=view.scene, radius=i.bfp_rad, length=0.1,
                                   rows=1, cols=100, cone_length=0.01,
                                   color=(1,1,0,1.0))
            m = MatrixTransform()
            m.translate((0, 0, 2*i.f))
            m.rotate(angle=angle, axis=axis)
            circ.transform = m

            # Plot polarizations
            if pol_dirs == None:
                pol_dirs = [i.bfp_pol_dir]
            for pol_dir in pol_dirs:
                for direction in [-1, 1]:
                    pol = visuals.MyArrow(parent=view.scene, length=i.bfp_rad/2, radius=0.25)
                    m = MatrixTransform()
                    m.rotate(angle=direction*90, axis=(0,1,0))
                    phi_angle = np.degrees(np.arctan2(pol_dir[1], pol_dir[0]))
                    m.rotate(angle=phi_angle, axis=(0,0,1))
                    m.translate((0, 0, 2*i.f))
                    m.rotate(angle=angle, axis=axis)
                    pol.transform = m
                    
        # Plot dipole
        dip = visuals.MyArrow(parent=view.scene, length=2, radius=0.25)
        m = MatrixTransform()
        m.rotate(angle=-45, axis=(1,1,0))
        dip.transform = m
        
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
