import numpy as np

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

def plot_sphere(filename, title, directions, data, display='save', random_colors=False, show_edges=False):
    from scipy.spatial import ConvexHull
    import vispy
    from vispy import scene
    from vispy.visuals.transforms import STTransform
    from vispy.geometry import MeshData
    
    # from vispy.gloo import Program
    vispy.use('glfw')

    canvas = scene.SceneCanvas(keys='interactive', bgcolor='white',
                               size=(1000, 1000), show=True)

    view = canvas.central_widget.add_view()
    view.camera = 'turntable'
    view.camera.azimuth = 135
    

    # Process Data
    # Clean this!
    def tp2xyz(tp):
        r = 0.8
        return [r*np.sin(tp[0])*np.cos(tp[1]),
                r*np.sin(tp[0])*np.sin(tp[1]),
                r*np.cos(tp[0])]
    
    points = np.apply_along_axis(tp2xyz, 1, directions)
    ch = ConvexHull(points)
    mesh = MeshData(vertices=ch.points, faces=ch.simplices)

    data_norm = data/(2*np.percentile(data, 90))
    def color_map(z):
        return (.75+.25*z[0],.25+.75*z[0],.25+.75*z[0],1.0)

    z = np.expand_dims(data_norm, 1)
    color = np.apply_along_axis(color_map, 1, z)
 
    # Random colors
    if random_colors:
        color = np.random.uniform(size=(mesh.get_vertices().shape[0], 4))
        color[:,3] = 1

    # Edges
    if show_edges:
        edge_color = 'black'
    else:
        edge_color = None
        
    MySphere(parent=view.scene, radius=0.8, vertex_colors=color,
                      edge_color=edge_color, mesh=mesh)
    MyXYZAxis(parent=view.scene, origin=[1.5,0,1.5], length=0.3)
    vispy.scene.Text(title, parent=view.scene, font_size=16, pos=(0,0,1.2))    

    
    # Display
    if(display=='save'):
        import scipy.misc
        scipy.misc.imsave(filename, canvas.render())
        
    elif(display=='interact'):
        canvas.app.run()

# MySphereVisual (move this to a separate file later)
from vispy.geometry import create_sphere
from vispy.visuals.mesh import MeshVisual
from vispy.visuals import CompoundVisual
from vispy.scene.visuals import create_visual_node

class MySphereVisual(CompoundVisual):
    
    def __init__(self, radius=1.0, directions=None,
                 edge_color='black', vertex_colors=None, mesh=None):
        
        self._mesh = MeshVisual(vertices=mesh.get_vertices(),
                                faces=mesh.get_faces(),
                                vertex_colors=vertex_colors)
                                
                                
        if edge_color:
            self._border = MeshVisual(vertices=mesh.get_vertices(),
                                      faces=mesh.get_edges(),
                                      color=edge_color, mode='lines')
        else:
            self._border = MeshVisual()

        CompoundVisual.__init__(self, [self._mesh, self._border])
        self.mesh.set_gl_state(polygon_offset_fill=True,
                               polygon_offset=(1, 1), depth_test=True)

    @property
    def mesh(self):
        """The vispy.visuals.MeshVisual that used to fil in.
        """
        return self._mesh

    @property
    def border(self):
        """The vispy.visuals.MeshVisual that used to draw the border.
        """
        return self._border

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
