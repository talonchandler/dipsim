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

def plot_sphere(filename, directions, data, display='save'):
    from glumpy import app, gl
    from glumpy.graphics.text import FontManager
    from glumpy.graphics.collections import TriangleCollection, PathCollection, GlyphCollection
    from glumpy.transforms import Position, Trackball, Viewport
    from glumpy.geometry import primitives
    from scipy.spatial import ConvexHull

    # Shaders
    vert = """
    //attribute vec4 a_color;
    varying vec4  v_color; 
    void main()
    {
        fetch_uniforms();
        v_color = a_color;
        gl_Position = <transform(position)>;
    }
    """
    frag = """
    varying vec4 v_color;
    void main(void)
    {
        gl_FragColor = v_color;
    }
    """
    
    # App Setup
    n = 2000
    app.use("glfw")
    window = app.Window(color=(1,1,1,1), width=n, height=n)
    
    @window.event
    def on_draw(dt):
        window.clear()
        cells.draw()
        #outlines.draw()
        axes.draw()        
        labels.draw()
        if(display=='save'):
            from glumpy.ext import png
            framebuffer = np.zeros((window.height, window.width * 3),
                                   dtype=np.uint8)
            gl.glReadPixels(0, 0, window.width, window.height,
                            gl.GL_RGB, gl.GL_UNSIGNED_BYTE, framebuffer)        
            png.from_array(framebuffer, 'RGB').save(filename)

    @window.event
    def on_init():
        transform.theta = 45
        transform.phi = 135
        transform.zoom = 16
        
        gl.glLineWidth(2)
        gl.glEnable(gl.GL_DEPTH_TEST)
        gl.glEnable(gl.GL_LINE_SMOOTH)

    @window.event
    def on_key_press(key, modifiers):
        if key == app.window.key.SPACE:
            on_init()

    transform = Trackball(Position())            
    viewport = Viewport(size=(1.0, 1.0))

    # Setup Collections
    cells = TriangleCollection("raw", user_dtype=[('a_color', (np.float32, 4), '!local', (0,0,0,1))], vertex=vert, fragment=frag, transform=transform)
    outlines = PathCollection("raw", transform=transform, color='shared')
    axes = PathCollection("raw", transform=transform, color='shared', viewport=viewport)    
    labels = GlyphCollection(transform=transform)

    # Process Data
    def tp2xyz(tp):
        return [np.sin(tp[0])*np.cos(tp[1]),
                np.sin(tp[0])*np.sin(tp[1]),
                np.cos(tp[0])]
    
    points = np.apply_along_axis(tp2xyz, 1, directions)
    de = ConvexHull(points)
    data_norm = data/(2*np.percentile(data, 90))
    
    for i, simplex in enumerate(de.simplices):
        def color_map(z):
            return (.75+.25*z[0],.25+.75*z[0],.25+.75*z[0],1.0)

        z = np.expand_dims(data_norm[simplex], 1)
        color = np.apply_along_axis(color_map, 1, z)
        
        V = de.points[simplex]
        I = np.zeros((len(V)-2,3))
        I[:,1] = 1 + np.arange(len(I))
        I[:,2] = 1 + I[:,1]

        cells.append(V, I.ravel(), a_color=color)        
        outlines.append(V, color=(0, 0, 0, 1), closed=True)

    # Axes
    X = np.array(([0.,0.,0.],[1.2,0.,0.]))
    Y = np.array(([0.,0.,0.],[0.,1.2,0]))
    Z = np.array(([0.,0.,0.],[0.,0.,1.2]))
    axes.append(X, color=(0,0,0,1), closed=False)
    axes.append(Y, color=(0,0,0,1), closed=False)
    axes.append(Z, color=(0,0,0,1), closed=False)    
    
    # Labels
    font = FontManager.get("OpenSans-Regular.ttf")
    scale = 0.002
    labels.append('x', font, scale=scale, origin=(1.3,0,0), direction=(-1,1,0))
    labels.append('y', font, scale=scale, origin=(0,1.3,0), direction=(-1,1,0))
    labels.append('z', font, scale=scale, origin=(0,0,1.3), direction=(-1,1,0))
    
    window.attach(outlines["transform"])
    window.attach(outlines["viewport"])
    window.attach(axes["viewport"])
    window.attach(labels["viewport"])    

    # Display
    if(display=='save'):
        app.run(framecount=1)
    elif(display=='interact'):
        app.run()
