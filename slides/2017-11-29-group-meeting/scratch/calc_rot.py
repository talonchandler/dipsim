import sympy.matrices.dense
from sympy import Symbol, simplify
theta = Symbol("theta")
phi = Symbol("phi")
Rytheta = sympy.dense.rot_axis2(-theta)
Rzphi = sympy.dense.rot_axis3(-phi)

print(Rytheta)
print(Rzphi)
print(simplify(Rzphi*Rytheta))


