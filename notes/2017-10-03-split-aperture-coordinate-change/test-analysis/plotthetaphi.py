import numpy as np
import util
import glob
files = glob.glob("./output/*.csv")

for file in files:
    print(file)
    data = np.genfromtxt(file, delimiter=',', skip_header=1,
                         names=['R-view x (um)', 'R-view y (um)', 'Theta [z polar angle] (rad)', 'Phi [x azimuth angle] (rad)'])

    # Create asymptote string
    asy_string = ''
    for entry in data:
        temp_str = 'dip_arrow(x_pos, y_pos, theta, phi, 10, (0,0,0));'
        temp_str = temp_str.replace('x_pos', str(entry['Rview_x_um']))
        temp_str = temp_str.replace('y_pos', str(entry['Rview_y_um']))
        temp_str = temp_str.replace('theta', str(entry['Theta_z_polar_angle_rad']))
        temp_str = temp_str.replace('phi', str(entry['Phi_x_azimuth_angle_rad']))        
        asy_string += temp_str
    
    # Generate asymptote plot
    util.draw_scene(asy_string, filename=file.replace('.csv', '.pdf'), save_file=True)
