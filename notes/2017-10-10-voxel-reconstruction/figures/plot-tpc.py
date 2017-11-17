import numpy as np
import h5py
import my_draw_scene

h5f = h5py.File('data.h5','r')
data = h5f['result'][:]
h5f.close()

# Append indices
nx, ny = data.shape[0], data.shape[1]
x = np.linspace(0, nx-1, nx)
y = np.linspace(0, ny-1, ny)
xv, yv = np.meshgrid(x, y)
xvv = np.expand_dims(xv, 2)
yvv = np.expand_dims(yv, 2)
data = np.concatenate([yvv, xvv, data], 2)

# Create asymptote string
def asy_string_wrapper(data):
    data = np.round(data, 3)
    temp_str = 'dip_arrow(x_pos, y_pos, theta, phi, 0.4*c, (0,0,0));'
    temp_str = temp_str.replace('x_pos', str(data[0]))
    temp_str = temp_str.replace('y_pos', str(data[1]))
    temp_str = temp_str.replace('theta', str(data[2]))
    temp_str = temp_str.replace('phi', str(data[3]))
    temp_str = temp_str.replace('c', str(data[4]))
    return temp_str

asy_string = ''
for datum in data.reshape((data.shape[0]*data.shape[1], data.shape[2])):
    asy_string += asy_string_wrapper(datum)


# Generate asymptote plot
my_draw_scene.draw_scene(asy_string, filename='output.pdf', save_file=True)
