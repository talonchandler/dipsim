#cython: boundscheck=False, wraparound=False, nonecheck=False
from libc.math cimport sin, cos, acos, atan2, exp, sqrt, M_PI
from libc.stdlib cimport malloc, free
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from cython_gsl cimport *

ctypedef double * double_ptr

cdef double calc_excitation_efficiency(double theta_lab, double phi_lab, double phi_pol, double theta_optical_axis) nogil:
    cdef double theta_opt, num, den, acos_arg
    theta_opt = acos(sin(theta_optical_axis)*cos(phi_lab)*sin(theta_lab) + cos(theta_optical_axis)*cos(theta_lab))
    num = cos(theta_optical_axis)*cos(phi_lab)*sin(theta_lab) - sin(theta_optical_axis)*cos(theta_lab)
    den = sqrt(1 - (sin(theta_optical_axis)*cos(phi_lab)*sin(theta_lab) + cos(theta_optical_axis)*cos(theta_lab))**2)

    # Prevent division by zero
    if den == 0:
        acos_arg = 1.0*num
    else:
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

cdef double calc_detection_efficiency(double theta_lab, double phi_lab, double theta_optical_axis, double alpha) nogil:
    cdef double theta_opt, A, B
    theta_opt = acos(sin(theta_optical_axis)*cos(phi_lab)*sin(theta_lab) + cos(theta_optical_axis)*cos(theta_lab))
    A = 0.25 - 0.375*cos(alpha) + 0.125*(cos(alpha)**3)
    B = 0.1875*(cos(alpha) - (cos(alpha)**3))
    return 2.0*(A + B*(sin(theta_opt)**2))

# Intensity calculation for a single fluorophore
cpdef double intensity(double theta_lab, double phi_lab, double phi_pol,
               double ill_theta_optical_axis, double det_theta_optical_axis,
               double det_alpha, double c) nogil:
    cdef double exc_eff, det_eff
    exc_eff = calc_excitation_efficiency(theta_lab, phi_lab, phi_pol, ill_theta_optical_axis)
    det_eff = calc_detection_efficiency(theta_lab, phi_lab, det_theta_optical_axis, det_alpha)
    return c*exc_eff*det_eff

# Normalization constant in Watson distributions
cdef double calc_norm(double kappa) nogil:
    return 1.0/(4*M_PI*gsl_sf_hyperg_1F1(0.5, 1.5, kappa))

# Fluorophore weight in Watson distribtuion
cdef double calc_weight(double theta1, double phi1, double theta2, double phi2, double kappa) nogil:
    cdef double arg, x1, y1, z1, x2, y2, z2
    x1 = sin(theta1)*cos(phi1)
    y1 = sin(theta1)*sin(phi1)
    z1 = cos(theta1)
    x2 = sin(theta2)*cos(phi2)
    y2 = sin(theta2)*sin(phi2)
    z2 = cos(theta2)
    arg = kappa*((x1*x2+y1*y2+z1*z2)**2)
    return exp(arg)

# Inner integrand for intensity from distribution of fluorophores
cdef double inner_integrand(double theta_lab, void *params) nogil:
    cdef double *cast_params = <double *> params;
    cdef double phi_lab = (<double_ptr> params)[0]
    cdef double phi_pol = (<double_ptr> params)[1]
    cdef double ill_theta_optical_axis = (<double_ptr> params)[2]
    cdef double det_theta_optical_axis = (<double_ptr> params)[3]
    cdef double det_alpha = (<double_ptr> params)[4]
    cdef double c = (<double_ptr> params)[5]
    cdef double theta_mean_flu = (<double_ptr> params)[6]
    cdef double phi_mean_flu = (<double_ptr> params)[7]
    cdef double kappa_flu = (<double_ptr> params)[8]
              
    cdef double I, norm, weight, jacobian
    I = intensity(theta_lab, phi_lab, phi_pol, ill_theta_optical_axis, det_theta_optical_axis, det_alpha, c)
    norm = calc_norm(kappa_flu)
    weight = calc_weight(theta_lab, phi_lab, theta_mean_flu, phi_mean_flu, kappa_flu)
    jacobian = sin(theta_lab)
    return I*norm*weight*jacobian

# Outer integrand for intensity from distribution of fluorophores
cdef double outer_integrand(double phi_lab, void *params) nogil:
    cdef double *cast_params = <double *> params;
    cast_params[0] = phi_lab;
    
    cdef gsl_integration_workspace *w
    cdef gsl_function F
    cdef double result, error
    w = gsl_integration_workspace_alloc(1000)
    F.function = &inner_integrand
    F.params = params

    gsl_integration_qag(&F, 0, M_PI, 0.01, 0, 1000, GSL_INTEG_GAUSS15, w, &result, &error)
    gsl_integration_workspace_free(w);
    return result

# Intensity integral for distribution of fluorophores
def intensity_from_dist(double phi_pol, double ill_theta_optical_axis,
                        double det_theta_optical_axis, double det_alpha, double c,
                        double theta_mean_flu, double phi_mean_flu, double kappa_flu):
    cdef double *params = <double *>PyMem_Malloc(10*sizeof(double))
    params[1] = phi_pol
    params[2] = ill_theta_optical_axis
    params[3] = det_theta_optical_axis
    params[4] = det_alpha
    params[5] = c
    params[6] = theta_mean_flu
    params[7] = phi_mean_flu
    params[8] = kappa_flu

    cdef gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000)
    cdef gsl_function F
    F.function = &outer_integrand
    F.params = params

    cdef double result, error
    gsl_integration_qag(&F, -M_PI, M_PI, 0.01, 0, 1000, GSL_INTEG_GAUSS15, w, &result, &error)
    gsl_integration_workspace_free(w)
    PyMem_Free(params)

    return result
