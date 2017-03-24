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
    def __init__(self, illuminator, detector):
        self.illuminator = illuminator
        self.detector = detector

    def calc_induced_dipoles(self, fluorophores):
        for f in fluorophores:
            f.mu_ind = f.mu_em*np.dot(f.mu_abs, self.illuminator.E_eff)

    def calc_total_intensity(self, fluorophores):
        self.calc_induced_dipoles(fluorophores)
        # TODO Green's tensor integrated over area
        # For now sum over entire volume
        I = 0
        for f in fluorophores:
            I += np.linalg.norm(f.mu_ind)**2
        
        return I

    def calc_total_intensity_from_single_fluorophore(self, args):
        theta = args[0]
        phi = args[1]

        flu_dir = np.array([np.sin(theta)*np.cos(phi),
                           np.sin(theta)*np.sin(phi),
                           np.cos(theta)])

        flu = fluorophore.Fluorophore(position=np.array([0, 0, 0]),
                                      mu_abs=flu_dir,
                                      mu_em=flu_dir)

        I = self.calc_total_intensity([flu])
        return I

    
    def plot_intensities_from_single_fluorophore(self, filename, n=50, **kwargs):
        directions = util.fibonacci_sphere(n)
        print('Generating data for microscope: '+filename)
        I = np.apply_along_axis(self.calc_total_intensity_from_single_fluorophore,
                                      1, directions)
        print('Plotting data for microscope: '+filename)
        util.plot_sphere(filename, directions=directions, data=I/I.max(), **kwargs)

    def draw_scene(self, filename, interact=False, my_ax=None, dpi=500, vis_px=1000):
        vispy.use('glfw')
        vis = vispy.scene.visuals
        
        # Setup viewing window
        canvas = vispy.scene.SceneCanvas(keys='interactive', bgcolor='white',
                                         size=(vis_px, vis_px), show=interact)
        cam = vispy.scene.cameras.turntable.TurntableCamera
        my_cam = cam(fov=0, azimuth=135, scale_factor=23,
                     center=(0, 0, self.illuminator.f))
        view = canvas.central_widget.add_view(camera=my_cam)

        # Plot dipole
        dip = visuals.MyArrow(parent=view.scene, length=2)
        m = MatrixTransform()
        m.rotate(angle=-45, axis=(1,1,0))
        dip.transform = m
        
        # Plot illuminator
        i = self.illuminator
        pol = visuals.MyArrow(parent=view.scene, length=1.5)
        m = MatrixTransform()
        m.rotate(angle=90, axis=(0,1,0))
        m.translate((0,0,2*i.f+0.1))
        pol.transform = m
        
        circ = vis.Ellipse(parent=view.scene, center=(0,0), radius=i.bfp_rad,
                           color='yellow', border_color='black', border_width=vis_px/150)
        circ.transform = STTransform(translate=(0, 0, 2*i.f))
        
        # Plot lens
        lens2 = visuals.MyLens(parent=view.scene)
        lens2.transform = STTransform(scale=(10, 10, 10),
                                      translate=(0, 0, i.f))
        lens_out = vis.Ellipse(parent=view.scene, center=(0,0), radius=0.5,
                               color=None, border_color='black', border_width=vis_px/150)
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
            f.savefig(filename, dpi=dpi)
            return ax

