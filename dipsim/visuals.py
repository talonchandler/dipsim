# Custom vispy classes

from vispy.geometry import create_sphere
from vispy.visuals.mesh import MeshVisual
from vispy.visuals import CompoundVisual
from vispy.scene.visuals import create_visual_node

class MySphereVisual(CompoundVisual):
    
    def __init__(self, radius=1.0, directions=None,
                 vertex_colors=None, mesh=None):
        
        self._mesh = MeshVisual(vertices=mesh.get_vertices(),
                                faces=mesh.get_faces(),
                                vertex_colors=vertex_colors)

        CompoundVisual.__init__(self, [self._mesh])
        self.mesh.set_gl_state(depth_test=True)

    @property
    def mesh(self):
        """The vispy.visuals.MeshVisual that used to fil in.
        """
        return self._mesh

MySphere = create_visual_node(MySphereVisual)

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
