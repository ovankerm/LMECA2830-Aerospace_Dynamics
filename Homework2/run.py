import numpy as np
from matplotlib import pyplot as plt

# ------ DATA ------
S = 2.01
b = 3.24
c_bar = S/b
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

t_star = c_bar/(2 * u)

Cy_b = -0.2833
Cl_b = -0.1031
Cn_b = 0.0344

Cy_p = 0
Cl_p = -0.218
Cn_p = -0.0328

Cy_r = 0
Cl_r = 0.0731
Cn_r = -0.0116


# ------ FUNCTIONS ------
def print_eigval(eigval, eigvect, adim_index = 3):
    real = False
    if eigval.imag == 0:
        eigval = eigval.real
        eigvect = eigvect.real
        real = True

    print(f"Eigenvalue {eigval} with eigenvector {eigvect}")

    if real:
        print("T_0.5 = %.4f"%(np.log(2)/np.abs(eigval)))
    else:
        print("T_0.5 = %.4f"%(np.log(2)/np.abs(eigval.real)))
        print("T = %.4f"%(2 * np.pi/np.abs(eigval.imag)))
        print("N_0.5 = %.4f"%(np.log(2) * np.abs(eigval.imag)/(np.abs(eigval.real) * 2 * np.pi)))

    psi = eigvect[2]/eigval
    y = (u * psi + eigvect[0])/eigval
    print(psi)
    var = np.array([eigvect[0]/u, eigvect[1]/(2 * u/b), eigvect[2]/(2 * u/b), eigvect[3], psi, y/(u * t_star)])
    var /= var[adim_index]
    print(f"beta = {np.abs(var[0])} {np.degrees(np.angle(var[0]))} \n p = {np.abs(var[1])} {np.degrees(np.angle(var[1]))} \n r = {np.abs(var[2])} {np.degrees(np.angle(var[2]))} \n phi = {np.abs(var[3])} {np.degrees(np.angle(var[3]))} \n psi = {np.abs(var[4])} {np.degrees(np.angle(var[4]))} \n y_0 = {np.abs(var[5])} {np.degrees(np.angle(var[5]))} \n")


print()

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

print("------ AERO DERIVATIVES ------")
print("Yv = %.4f Lv = %.4f Nv = %.4f"%(Y_v, L_v, N_v))
print("Yp = %.4f Lp = %.4f Np = %.4f"%(Y_p, L_p, N_p))
print("Yr = %.4f Lr = %.4f Nr = %.4f"%(Y_r, L_r, N_r))
print()


# ------ MATRIX ------

print(N_v)

A = np.array([[Y_v/m_gross, Y_p/m_gross, Y_r/m_gross - u, g],
              [I_zz_p * L_v - I_xz_p * N_v, I_zz_p * L_p - I_xz_p * N_p, I_zz_p * L_r - I_xz_p * N_r, 0],
              [I_xx_p * N_v - I_xz_p * L_v, I_xx_p * N_p - I_xz_p * L_p, I_xx_p * N_r - I_xz_p * L_r, 0],
              [0, 1, 0, 0]])

print("------ A MATRIX ------")
print(A)
print()


# ------ POLYNOMIAL ------
eigvals, eigvect = np.linalg.eig(A)
p = np.polynomial.Polynomial(np.real(np.polynomial.polynomial.polyfromroots(eigvals)))
print("------ CHARACTERISTIC POLYNOMIAL ------")
print(p)
print()

print("------ EIGENVALUES ------")
# print(eigvals)
# print()

# print("------ EIGENVECTORS ------")
# for i in range(len(eigvals)):
#     print(f"Eigval : {eigvals[i]} with eigvect : {eigvect[:, i]}")
# print()


indices = [3, 3, 3, 4]
for i in range(4) : print_eigval(eigvals[i], eigvect[:, i], indices[i])



