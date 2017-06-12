import numpy as np

alpha = np.pi/2

A = (1.0/6.0) - 0.25*np.cos(alpha) + (1.0/12.0)*(np.cos(alpha)**3)
B = 0.125*np.cos(alpha) - 0.125*(np.cos(alpha)**3)

print(A)
print(B)
print(2*A + 2*B)
