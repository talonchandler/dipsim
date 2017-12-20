import numpy as np
from mpmath import *
from sympy import *
from sympy.matrices.dense import *
import sympy.functions.special.spherical_harmonics as sh
import functools

print("Working...")

# Symbols
theta = Symbol('theta', real=True)
phi = Symbol('phi', real=True)
Theta = Symbol('Theta', real=True)
Phi = Symbol('Phi', real=True)

phi_pol = Symbol('phi_pol', real=True)
A = Symbol('A')
B = Symbol('B')
C = Symbol('C')
D = Symbol('D')

# Takes an expression and returns list of coefficients in (l, m, coeff) tuples.
def ylm_coeffs(model, max_l=4):
    coeffs = []
    for n in range(0, max_l+2, 2):
        coeff_row = []        
        for m in range(-n, n+1):
            print("Integrating: "+ str(n) + ', ' + str(m))            
            ynm = conjugate(sh.Ynm(n, m, theta, phi).expand(func=True))
            theta_int = integrate(expand(sin(theta)*ynm*model), (theta, 0, pi)) # theta integral
            final_int = simplify(integrate(expand_trig(theta_int), (phi, 0, 2*pi))) # phi integral
            coeff_row.append(final_int)
        coeffs.append(coeff_row)
    return np.array(coeffs)

# Other models
#model1 = 1 - sin(theta)**2*cos(phi)**2 # x detector
#model1 = 1 - sin(theta)**2*sin(phi)**2 # y detector
#model1 = sin(theta)**2 # z detector
#model1 = 1 - (sin(Theta)*cos(Phi)*sin(theta)*cos(phi) + sin(Theta)*sin(Phi)*sin(theta)*sin(phi) + cos(Theta)*cos(theta))**2 # Theta, Phi detector
#model2 = 1
#model =  D*(A + B*(sin(theta)**2) + C*(sin(theta)**2)*(cos(phi - phi_pol)**2))*(A + B*(sin(theta)**2)) # epipol-epi

# diSPIM model
model1 = (sin(theta)**2)*(cos(phi - phi_pol)**2)*2*(A + B*(cos(theta)**2 + sin(theta)**2*sin(phi)**2)) # orthopol-epi (illumination from z - detection on x)
model2 = ((sin(theta)*sin(phi)*sin(phi_pol) - cos(theta)*cos(phi_pol))**2)*2*(A + B*(sin(theta)**2)) # orthopol-epi (illumination from x - detection on z)

c1 = ylm_coeffs(model1, max_l=6)
c2 = ylm_coeffs(model2, max_l=6)
print("ill z, det x, phi_pol (x -> y)")
print(c1)
print("ill x, det , phi_pol (-z -> y)")
print(c2)
