#cython: boundscheck=False, wraparound=False, nonecheck=False
from libc.math cimport sin, cos, acos, sqrt, M_PI

def cy_calc_collection_efficiency(double theta_lab, double phi_lab, double theta_optical_axis, double alpha):
    cdef double theta_opt, A, B
    theta_opt = acos(sin(theta_optical_axis)*cos(phi_lab)*sin(theta_lab) + cos(theta_optical_axis)*cos(theta_lab))
    A = 0.25 - 0.375*cos(alpha) + 0.125*(cos(alpha)**3)
    B = 0.1875*(cos(alpha) - (cos(alpha)**3))
    return 2.0*(A + B*(sin(theta_opt)**2))

def cy_calc_excitation_efficiency(double theta_lab, double phi_lab, double phi_pol, double theta_optical_axis):
    cdef double theta_opt, num, den, acos_arg
    theta_opt = acos(sin(theta_optical_axis)*cos(phi_lab)*sin(theta_lab) + cos(theta_optical_axis)*cos(theta_lab))
    num = cos(theta_optical_axis)*cos(phi_lab)*sin(theta_lab) - sin(theta_optical_axis)*cos(theta_lab)
    den = sqrt(1 - (sin(theta_optical_axis)*cos(phi_lab)*sin(theta_lab) + cos(theta_optical_axis)*cos(theta_lab))**2)
    acos_arg = num/den

    # Prevent nan from floating point errors
    if acos_arg > 1.0:
        acos_arg = 1.0
    if acos_arg < -1.0:
        acos_arg = -1.0

    # Check sign
    if phi_lab < M_PI and phi_lab > 0:
        phi_opt = acos(acos_arg)
    else:
        phi_opt = -acos(acos_arg)
        
    return (sin(theta_opt)*cos(phi_opt - phi_pol))**2
