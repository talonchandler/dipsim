# Custom vispy classes
from dipsim import util
import numpy as np
from vispy.geometry import create_sphere, MeshData
from vispy.visuals.mesh import MeshVisual
from vispy.visuals import CompoundVisual
from vispy.scene.visuals import create_visual_node
import matplotlib
import matplotlib.pyplot as plt
from vispy.visuals.transforms import (STTransform, LogTransform,
                                      MatrixTransform, PolarTransform)


class MySphereVisual(CompoundVisual):
    
    def __init__(self, radius=1.0, directions=None, colors=None):
                 
        # Convert spherical to cartesian
        points = np.apply_along_axis(util.tp2xyz, 1, directions)

        # Create mesh
        import scipy.spatial
        ch = scipy.spatial.ConvexHull(points)
        mesh = MeshData(vertices=ch.points, faces=ch.simplices)

        self._mesh = MeshVisual(vertices=mesh.get_vertices(),
                                faces=mesh.get_faces(),
                                vertex_colors=colors)

        CompoundVisual.__init__(self, [self._mesh])
        self.mesh.set_gl_state(depth_test=True)

    @property
    def mesh(self):
        """The vispy.visuals.MeshVisual that used to fil in.
        """
        return self._mesh

MySphere = create_visual_node(MySphereVisual)

class MyLensVisual(CompoundVisual):
    
    def __init__(self, radius=1, max_theta=np.pi/6, n_theta=16, n_phi=64, color=(0.5,0.5,0.5,1.0)):
        # Create lens points based on inputs
        phi = np.linspace(0, 2*np.pi, n_phi, endpoint=False)
        theta = np.linspace(max_theta, 0, n_theta, endpoint=False)        
        directions = np.array(np.meshgrid(theta, phi)).T.reshape(-1, 2)
        t_points = np.apply_along_axis(util.tp2xyz, 1, directions)
        t_points[:,2] -= np.cos(max_theta)
        b_points = np.copy(t_points)
        b_points[:,2] = -t_points[:,2]
        points = np.concatenate((t_points, b_points))
        
        # Create mesh
        import scipy.spatial
        ch = scipy.spatial.ConvexHull(points)
        mesh = MeshData(vertices=ch.points, faces=ch.simplices)
        
        self._mesh = MeshVisual(vertices=mesh.get_vertices(),
                                faces=mesh.get_faces(),
                                color=color)

        CompoundVisual.__init__(self, [self._mesh])
        self.mesh.set_gl_state(depth_test=True)

    @property
    def mesh(self):
        """The vispy.visuals.MeshVisual that used to fil in.
        """
        return self._mesh

MyLens = create_visual_node(MyLensVisual)


class MyArrowVisual(CompoundVisual):
    
    def __init__(self, rows=30, cols=30, radius=0.1, length=1, color='black', cone_length=0.5, cone_radius=None):

        import vispy.geometry.generation as gen
        mesh = gen.create_arrow(rows=rows, cols=cols, radius=radius, length=length, cone_length=cone_length, cone_radius=cone_radius)
        
        self._mesh = MeshVisual(vertices=mesh.get_vertices(),
                                faces=mesh.get_faces(),
                                color=color)

        CompoundVisual.__init__(self, [self._mesh])
        self.mesh.set_gl_state(depth_test=False)

    @property
    def mesh(self):
        """The vispy.visuals.MeshVisual that used to fil in.
        """
        return self._mesh

MyArrow = create_visual_node(MyArrowVisual)

class MyLineVisual(CompoundVisual):
    
    def __init__(self, rows=10, cols=30, radius=0.1, length=1, color='black'):

        import vispy.geometry.generation as gen
        mesh = gen.create_cylinder(rows=rows, cols=cols, radius=[radius, radius], length=length)
        
        self._mesh = MeshVisual(vertices=mesh.get_vertices(),
                                faces=mesh.get_faces(),
                                color=color)

        CompoundVisual.__init__(self, [self._mesh])
        self.mesh.set_gl_state(depth_test=True)

    @property
    def mesh(self):
        """The vispy.visuals.MeshVisual that used to fil in.
        """
        return self._mesh

MyLine = create_visual_node(MyLineVisual)

from vispy.visuals.line import LineVisual
from vispy.visuals.text import TextVisual
class MyXYZAxisVisual(CompoundVisual):
    """
    Simple 3D axis for indicating coordinate system orientation. Axes are
    x=red, y=green, z=blue.
    """
    def __init__(self, origin=[0,0,0], length=1):
        verts = origin + np.array([[0, 0, 0],
                                   [length, 0, 0],
                                   [0, 0, 0],
                                   [0, length, 0],
                                   [0, 0, 0],
                                   [0, 0, length]])

        line = LineVisual(pos=verts, color=np.array([0, 0, 0, 1]),
                          connect='segments', method='gl')

        x = TextVisual('x', font_size=12, pos=origin + np.array([1.25*length,0,0]))
        y = TextVisual('y', font_size=12, pos=origin + np.array([0,1.25*length,0]))
        z = TextVisual('z', font_size=12, pos=origin + np.array([0,0,1.25*length]))

        CompoundVisual.__init__(self, [line, x, y, z])

MyXYZAxis = create_visual_node(MyXYZAxisVisual)        
