import numpy as np
import util
import matplotlib.pyplot as plt
import os

def fibonacci_sphere(n):
    # Returns "equally" spaced points on a unit sphere in spherical coordinates.
    # http://stackoverflow.com/a/26127012/5854689
    z = np.linspace(1 - 1/n, -1 + 1/n, num=n) 
    theta = np.arccos(z)
    phi = np.mod((np.pi*(3.0 - np.sqrt(5.0)))*np.arange(n), 2*np.pi) - np.pi
    return np.vstack((theta, phi)).T

def phi_prime(theta, phi, alpha=np.pi/6):
    # Returns phi coordinate in a frame rotated by a right-handed rotation by
    # angle alpha about the positive y-axis. 
    if theta == 0:
        return 0
    num = np.cos(alpha)*np.cos(phi)*np.sin(theta) - np.sin(alpha)*np.cos(theta)
    den = np.sqrt(1 - (np.sin(alpha)*np.cos(phi)*np.sin(theta) + np.cos(alpha)*np.cos(theta))**2)
    if phi < np.pi and phi > 0:
        return np.arccos(num/den)
    elif phi < 0 and phi > -np.pi:
        return -np.arccos(num/den)
    else:
        return 0

def thetaphi2phiab(theta, phi, alpha=np.pi/6):
    # Converts theta, phi coordinates into phi_a, phi_b coordinates
    phi_a = phi_prime(theta, phi, alpha=alpha)
    phi_b = phi_prime(theta, phi, alpha=-alpha)
    return np.array([phi_a, phi_b])

def generate_LUTs(n=1e4, alpha=np.pi/6):
    # Generates LUTs in both coordinate systems
    thetaphiLUT = fibonacci_sphere(n)
    phiabLUT = []
    for direction in thetaphiLUT:
        phiabLUT.append(thetaphi2phiab(*direction, alpha=alpha))
    phiabLUT = np.array(phiabLUT)
    return thetaphiLUT, phiabLUT

def phiab2thetaphi(phia, phib, thetaphiLUT=None, phiabLUT=None):
    # Converts phi_a, phi_b coordinates into theta, phi using lookup tables.
    if phiabLUT is None:
        thetaphiLUT, phiabLUT = generate_LUTs()
    target = np.array([phia, phib])
    min_ind = np.argmin(np.linalg.norm(target - phiabLUT, axis=1))
    return thetaphiLUT[min_ind]

# Find .xlsx in current directory
import glob
files = glob.glob("./input/*.xlsx")

# Open raw data
import xlrd
book = xlrd.open_workbook(files[0])

# Generate LUT (time consuming)
print('Computing LUT.')
thetaphiLUT, phiabLUT = generate_LUTs(n=5e5)

# Organize data
for sheet_name in book.sheet_names():
    print('Processing: ', sheet_name)
    sheet = book.sheet_by_name(sheet_name)
    r_view_x = sheet.col_values(0, start_rowx=2)    
    r_view_y = sheet.col_values(1, start_rowx=2)
    r_view_ori = sheet.col_values(2, start_rowx=2)
    l_view_x = sheet.col_values(5, start_rowx=2)
    l_view_y = sheet.col_values(6, start_rowx=2)        
    l_view_ori = sheet.col_values(7, start_rowx=2)
    theta = []
    phi = []
    
    for phi_l, phi_r in zip(l_view_ori, r_view_ori):
        # Compute new coordinates
        tp = phiab2thetaphi(phi_l, phi_r, thetaphiLUT=thetaphiLUT, phiabLUT=phiabLUT)
        theta.append(tp[0])
        phi.append(tp[1])

    # Save data to csv
    save_data = np.vstack([r_view_x, r_view_y, r_view_ori, l_view_x, l_view_y, l_view_ori, theta, phi]).T
    header = 'R-view x (um), R-view y (um), R-view ori (rad), L-view x (um), L-view y (um), L-view ori (rad), Theta [z polar angle] (rad), Phi [x azimuth angle] (rad)'
    np.savetxt('./output/'+sheet_name.replace(' ', '')+'.csv', save_data, comments='', delimiter=',', header=header)

os.system('say "done"')
