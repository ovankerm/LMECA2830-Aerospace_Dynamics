import numpy as np
from matplotlib import pyplot as plt

# ------ QUESTION 1 -----
print("------ Answers for question 1 ------")
# ----PARAMS----
b = 13.62                  # [m]
S = 46.45                  # [m^2]
C_L_max = 1.8              # [-]
m = 13865                  # [kg]
n_max = 8                  # [-]
TSFC = 24.6e-6             # [kg/(N s)]
rho = 0.8189               # [kg/m^3]
g = 9.81                   # [m/s^2]
e = 0.74                   # [-]
C_D_0 = 0.026598           # [-]
delta_x = 600e3            # [km]

AR = b*b/S                 # [-]
k = 1./(np.pi * e * AR)     # [-]

#----POLAR----
C_L = np.linspace(-0.5, C_L_max)
C_D = C_D_0 + k * np.power(C_L, 2)

fig, ax = plt.subplots(1, 1, figsize=(10, 9))
ax.plot(C_D, C_L, 'k', label='Polar')
ax.grid()
ax.set_title("F/A-18 Super Hornet (ideal) polar curve")
ax.set_xlabel(r"$C_D$", size=12)
ax.set_ylabel(r"$C_L$", size=12)
ax.scatter(C_D_0, 0, c='r', label=r"Zero lift drag $C_{D,0}$")
ax.legend()
# fig.savefig("Homework1/images/Hornet_polar.eps", format='eps')

#----CL opt----
C_L_opt = np.sqrt(C_D_0/(3 * k))
C_D_opt = C_D_0 + k * np.power(C_L_opt, 2)
print(f"Optimal (C_L, C_D) = ({C_L_opt}, {C_D_opt})")

W_end = g * m
W_init = np.power(np.power(g * m, 0.5) + g * TSFC/2 * C_D_opt * np.power(rho * S/(2 * C_L_opt), 0.5) * delta_x, 2)
m_init = W_init/g
m_fuel = m_init - m
print(f"Fuel needed : {m_fuel} kg")

V_init = np.power(2 * W_init/(rho * C_L_opt * S), 0.5)
V_end = np.power(2 * W_end/(rho * C_L_opt * S), 0.5)

print(f"Initial velocity : {V_init} m/s = {V_init * 3.6} km/h")
print(f"Final velocity : {V_end} m/s = {V_end * 3.6} km/h")

#----LOAD FACTOR----
W = (m + 4000)*g
V_star = np.sqrt(2 * n_max * W/(S * rho * C_L_max))

print(f'Critical velocity : {V_star} m/s = {3.6 * V_star} km/h')

V_1 = np.linspace(0, V_star)
V_2 = np.linspace(V_star, 1.5 * V_star)

n_1 = 0.5 * rho * V_1**2 * C_L_max * S/W
n_2 = n_max * np.ones_like(V_2)

fig, ax = plt.subplots(1, 1, figsize=(10, 9))
ax.plot(V_1, n_1, 'k', label='Stall limit')
ax.plot(V_2, n_2, 'k--', label='Structural limit')
ax.grid()
ax.set_title("Flight envelope")
ax.set_xlabel("V")
ax.set_ylabel("n")
ax.scatter(V_star,n_max, c='r', label='Critical velocity')
ax.legend()
# fig.savefig('Homework1/images/flight_env.eps', format='eps')


#----TURNING RADIUS----
R_min_1 = 2 * n_max * W/(S * rho * C_L_max * (n_max + 1) * g)

print(f"R_min_1 = {R_min_1}")

R_min_2 = 2 * n_max * W/(S * rho * C_L_max * (n_max - 1) * g)

print(f"R_min_2 = {R_min_2}")


print()

# ------ QUESTION 2 -----
print("------ Answers for question 2 ------")
#----PARAMS----
m_gross = 3800
h = 0
S = 22.76
b = 11.25
c_bar = S/b
h_nw = 0
CL_w_a = 4.62
C_mac_w = -0.053
alpha_no_lift = -0.2
e = 0.83
C_D_0 = 0.0229
S_t = 2.9
b_t = 3.25
c_bar_t = S_t/b_t
h_nt = 2.74
CL_t_a = 4.06
C_mac_t = 0
i_t = 0
eps_0 = 0.035797
deps_da = 0.48833
CL_t_deps = 2.436

AR = b*b/S                 
k = 1./(np.pi * e * AR)     

#----STATIC MARGIN----
V_H_old = (h_nt + h - h_nw) * S_t/S
a_old = CL_w_a + S_t/S * CL_t_a * (1 - deps_da)
h_n_old = h_nw + 1/a_old * (V_H_old * CL_t_a * (1 - deps_da))

print(f"static margin : {h_n_old}")

#----NEW STATIC MARGIN----
V_H_new = (h_nt + h - h_nw) * (0.8 * S_t)/S
a_new = CL_w_a + (0.8 * S_t)/S * CL_t_a * (1 - deps_da)
h_n_new = h_nw + 1/a_new * (V_H_new * CL_t_a * (1 - deps_da))

print(f"new static margin : {h_n_new}")

#----K CONSTANT----
K = ((h_n_old - h_nw) * a_new - V_H_new * CL_t_a * (1 - deps_da))/(CL_t_deps * (V_H_new - (0.8 * S_t)/S * (h_n_old - h_nw)))

print(f"K = {K}")

#----C_trim plot----
alpha_trim = np.linspace(-np.radians(5), np.radians(10), num=2)
C_m_inf = C_mac_w + V_H_old * CL_t_a * (i_t + eps_0) * (1 - (CL_t_a * S_t)/(a_old * S) * (1 - deps_da))
C_L_trim = (a_old * alpha_trim * ((h_n_old - h_nw) * S_t/S - V_H_old) - S_t/S * C_m_inf)/((h - h_nw) * S_t/S - V_H_old)

fig, ax = plt.subplots(1, 1, figsize=(10, 9))
ax.grid()
ax.plot(np.degrees(alpha_trim), C_L_trim, 'k-', label = 'trimmed')
ax.plot(np.degrees(alpha_trim), alpha_trim * a_old, 'k--', label = 'untrimmed')
ax.legend()
ax.set_xlabel(r'$\alpha_{trim}$ [°]', size=12)
ax.set_ylabel(r'$C_L$', size=12)
fig.suptitle('Trimmed curve')
# fig.savefig('Homework1/images/trimmed_curve.eps', format='eps')

#----MAX RANGE----
denom = a_old * CL_t_deps * ((h_n_old - h_nw) * S_t/S - V_H_old)
C_L_opt = np.sqrt(C_D_0/(3 * k))
V = np.sqrt(2 * g * m_gross/(rho * C_L_opt * S))
alpha_opt = np.degrees((CL_t_deps * (S_t/S * (h - h_nw) - V_H_old) * C_L_opt + S_t/S * CL_t_deps * C_m_inf)/denom)
delta_e_opt = np.degrees(-a_old * (C_m_inf + (h - h_n_old) * C_L_opt)/denom)

print(f'Optimal AoA = {alpha_opt}')
print(f'Optimal EAS : {V} m/s = {3.6 * V} km/h')
print(f'Optimal aileron deflection : {delta_e_opt}')


plt.show()