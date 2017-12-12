import numpy as np
from mpmath import *
from sympy import *
from sympy.matrices.dense import *
import sympy.functions.special.spherical_harmonics as sh
import functools

print("Working...")

# Symbols
theta = Symbol('theta')
phi = Symbol('phi')
phi_pol = Symbol('phi_pol')
A = Symbol('A')
B = Symbol('B')
C = Symbol('C')
D = Symbol('D')

# Takes an expression and returns list of coefficients in (l, m, coeff) tuples.
def ylm_coeffs(model, max_l=4):
    coeffs = []
    for n in range(0, max_l+2, 2):
        complex_coeffs = []
        for m in range(-n, n+1):
            print("Integrating: "+ str(n) + ', ' + str(m))            
            ynm = simplify(expand_trig(sh.Ynm(n, m, theta, phi).expand(func=True).rewrite(cos)))
            theta_int = integrate(expand(sin(theta)*ynm*model), (theta, 0, pi)) # theta integral
            final_int = simplify(integrate(expand_trig(theta_int), (phi, 0, 2*pi))) # phi integral
            complex_coeffs.append(final_int)
        coeffs.append(c2r_coeffs(complex_coeffs))
    return np.array(coeffs)

# Convert coefficients of complex Ylm to coefficients of real ylm
def c2r_coeffs(coeffs):
    real_coeffs = []
    mid_ind = floor(len(coeffs)/2)
    for i, coeff in enumerate(coeffs):
        if i < mid_ind:
            real_coeffs.append((I/sqrt(2))*(coeffs[i] + coeffs[-i-1]))
        elif i == mid_ind:
            real_coeffs.append(coeffs[i])
        elif i > mid_ind:
            real_coeffs.append((1/sqrt(2))*(-coeffs[i] + coeffs[-i-1]))
    return coeffs

#model = cos(theta)**2 + sin(theta)**2*sin(phi)**2 # x detector
model = cos(theta)**2 + sin(theta)**2*cos(phi)**2 # y detector
#model = sin(theta)**2 # z detector

#model =  D*(A + B*(sin(theta)**2) + C*(sin(theta)**2)*(cos(phi - phi_pol)**2))*(A + B*(sin(theta)**2)) # epipol-epi
#model = 2*(sin(theta)**2)*(cos(phi - phi_pol)**2)*(A + B*(1 - (cos(phi)**2)*(sin(theta)**2))) # orthopol-epi (illumination from z - detection on x)
#model = 2*(sin(theta)**2)*(cos(phi - phi_pol)**2)*(A + B*(1 - (cos(phi)**2)*(sin(theta)**2))) # TODO orthopol-epi (illumination from x - detection on z)
c = ylm_coeffs(model, max_l=6)
print(c)
