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

def plot_sphere(filename, theta, phi, data):
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np

    # The Cartesian coordinates of the unit sphere
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)

    # Normalize to [0, 1]
    fmin = np.percentile(data, 0) 
    fmax = np.percentile(data, 95) # Use the first histogram bin as the max
    data = (data - fmin)/(fmax - fmin)
    
    # Set the aspect ratio to 1 so our sphere looks spherical
    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot(111, projection='3d', aspect='equal')
    cax = ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=cm.seismic(data))
    
    m = cm.ScalarMappable(cmap=cm.seismic)
    m.set_array(data)
    plt.colorbar(m)
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')    

    # Turn off the axis planes
    ax.set_axis_off()
    fig.savefig(filename)
