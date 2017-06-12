import numpy as np
from mpmath import *
from sympy import *
from sympy.matrices.dense import *
import functools
import util

print("Working...")

# Symbols
Theta = Symbol('Theta')
Phi = Symbol('Phi')
theta = Symbol('theta')
phi = Symbol('phi')
n0 = Symbol('n0')
n1 = Symbol('n1')
phi_pol = Symbol('phi_pol')
alpha = Symbol('alpha')

# Set up integral
mu = Matrix([sin(Theta)*cos(Phi), sin(Theta)*sin(Phi), cos(Theta)]) # Dipole 
r = Matrix([sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]) # Position
R = rot_axis3(-phi)*rot_axis2(theta)*rot_axis3(phi) # Rotation matrix
apod = 1 # sqrt(n1)/(sqrt(n0*cos(theta))) # Energy conservation
P = Matrix([cos(phi_pol), sin(phi_pol), 0]) # Polarizer
J = sin(theta) # Jacobian
integrand = J*(P.dot(apod*R*r.cross(mu).cross(r))**2) # Integrand
tot_integrand = J*(r.cross(mu).cross(r)).dot((r.cross(mu).cross(r))) 

# Calculate I integral
int1 = integrate(expand(integrand), (phi, 0, 2*pi)) # phi integral
I = integrate(expand(int1), (theta, 0, alpha)) # theta integral

# Calculate I_tot integral
int1 = integrate(expand(tot_integrand), (phi, 0, 2*pi)) # phi integral
I_tot = integrate(expand(int1), (theta, 0, pi)) # theta integral

# Calculate fraction of total power detected
I_frac = I/I_tot

# Print expressions for each polarization in Fourkas' form
pols = [0, pi/4, pi/2, 3*pi/4]
for pol in pols:
    I_pol = trigsimp(I_frac.subs([(phi_pol, pol)]))
    B_term = sin(Theta)**2
    if pol == 0 or pol == pi/2:
        sub_in = sin(Phi)**2
        sub_out = (1 - cos(2*Phi))/2
        C_term = sin(Theta)**2*cos(2*Phi)
    else:
        sub_in = sin(Phi)*cos(Phi)
        sub_out = sin(2*Phi)/2
        C_term = sin(Theta)**2*sin(2*Phi)
    if pol == 0 or pol == pi/4:
        C_sign = 1
    else:
        C_sign = -1
    
    I_pol_sin2cos = collect(I_pol, sub_in).subs([(sub_in, sub_out)])
    I_pol_collect_C_term = collect(expand(I_pol_sin2cos), C_sign*C_term)
    C_sin = C_sign*I_pol_collect_C_term.args[-1].args[1] # Extract C coefficient
    C = simplify(C_sin.subs([(sin(alpha)**2, 1 - cos(alpha)**2)])) # Write in cos
    A_and_B_terms = collect(Add(*I_pol_collect_C_term.args[0:-1]), B_term)
    B = A_and_B_terms.args[-1].args[1] # Extract B coefficient
    A = Add(*A_and_B_terms.args[0:-1]) # Extract A coefficient

    # Express coefficients in of powers of cos(alpha)
    coeffs = [expand_trig(x.subs([(sin(alpha)**2, 1 - cos(alpha)**2)])) for x in [A, B, C]]
    
    my_I_pol = coeffs[0] + coeffs[1]*B_term + coeffs[2]*C_sign*C_term

    # Ensure that simplification didn't change the expression
    assert(simplify(expand(my_I_pol) - expand(I_pol)) == 0) 

    # Print the expression in plain text
    # print(my_I_pol)

    # Print the expression in latex
    init_printing(order='rev-lex')
    names = ['A', 'B', 'C']
    int_string = 'I_{'+str(int(np.rad2deg(float(pol))))+'}(\\Theta, \\Phi, \\alpha) =& I_{\\text{tot}}(t, t+\\tau)(A + B' + latex(B_term) + ' + C(' + latex(C_sign*C_term) + '))\\\\'
    print(int_string.replace('\\left (', '').replace(' \\right )', ''))
    coeff_strings = [name + ' =& ' + latex(coeff).replace('\\\\', '\\') + '\\\\' for name, coeff in zip(names, coeffs)]
    for s in coeff_strings:
        print(s.replace('\\left (', '').replace('\\right )', ''))
    print('\n')
