import numpy as np
from matplotlib import pyplot as plt

# ------ DATA ------
S = 2.01
b = 3.24
m_gross = 12.5
I_xx = 5.742
I_yy = 5.98977
I_zz = 11.45476
I_xz = 0.07243

I_xx_p = I_xx/(I_xx * I_zz - I_xz*I_xz)
I_zz_p = I_zz/(I_xx * I_zz - I_xz*I_xz)
I_xz_p = I_xz/(I_xx * I_zz - I_xz*I_xz)

M = 0.06
rho = 1.225
p = 101325
T = 288.15
g = 9.81

R = p/(rho * T)
gamma = 1.4

u = M * np.sqrt(gamma * R * T)

Cy_b = -0.2833
Cl_b = -0.1031
Cn_b = -0.0344

Cy_p = 0
Cl_p = -0.218
Cn_p = -0.0328

Cy_r = 0
Cl_r = 0.0731
Cn_r = -0.0116


# ------ AERODYNAMIC DERIVATIVES ------
factor = 0.5 * rho * u * S

Y_v = factor * Cy_b
L_v = b * factor * Cl_b
N_v = b * factor * Cn_b

Y_p = 0.5 * b * factor * Cy_p
L_p = 0.5 * b*b * factor * Cl_p
N_p = 0.5 * b*b * factor * Cn_p

Y_r = 0.5 * b * factor * Cy_r
L_r = 0.5 * b*b * factor * Cl_r
N_r = 0.5 * b*b * factor * Cn_r

# ------ MATRIX ------
A = np.array([[Y_v/m_gross, Y_p/m_gross, Y_r/m_gross - u, g],
              [I_zz_p * L_v - I_xz_p * N_v, I_zz_p * L_p - I_xz_p * N_p, I_zz_p * L_r - I_xz_p * N_r, 0],
              [I_xx_p * L_v - I_xz_p * N_v, I_xx_p * L_p - I_xz_p * N_p, I_xx_p * L_r - I_xz_p * N_r, 0],
              [0, 1, 0, 0]])


print(A)

print(np.linalg.eigvals(A))
print(np.polynomial.polynomial.polyfromroots(np.linalg.eigvals(A)))