from scipy.integrate import ode
import matplotlib.pyplot as plt
from alu_model import *
from models import *
from parameters import *
import numpy as np


def int_to_bool_list(num):
    return [bool(num & (1 << n)) for n in range(4)]


def list_to_int(xs):
    v = 0
    f = 1
    for x in xs:
        v += x*f
        f *= 2
    return v

rho_x = 0

params = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X

# Vhodne logicne vrednosti
M = ['left', 'right']
A0 = [0,1]
A1 = [0,1]

# simulation parameters
t_end = 1000
N = t_end


# sums = []
# x_axis = []
# y_axis = []
for m in M:
    # left shift
    m0 = 0
    m1 = 1

    # overwrite if right shift
    if m == 'right':
        m0 = 1
        m1 = 0

    print(f'MODE: {m} - {m1}{m0}')

    for a1 in A1: 
        for a0 in A0:
            """
            SHIFTER

            A0, A1, M0, M1, I0, I1, \
            L_A0, L_A1, L_M0, L_M1, L_I0, L_I1, \
            N_A0, N_A1, N_M0, N_M1, N_I0, N_I1, \
            S0, S1
            """
            
            Y0 = np.zeros(20)

            Y0[:4] = a0, a1, m0, m1
            # N_X nastavimo vse na 1 - st celic
            Y0[12:-2] = 1

            T = np.linspace(0, t_end, N)

            t1 = t_end
            dt = t_end/N
            T = np.arange(0,t1+dt,dt)
            Y = np.zeros([1+N,len(Y0)])
            Y[0,:] = Y0

            r = ode(shifter_model_ODE).set_integrator('zvode', method='bdf')
            r.set_initial_value(Y0, T[0]).set_f_params(params)
            
            i = 1
            while r.successful() and r.t < t1:
                Y[i,:] = r.integrate(r.t+dt)
                i += 1

            s0 = 1 if Y[:,-2][-1] > 2 else 0
            s1 = 1 if Y[:,-1][-1] > 2 else 0

            # sums[0].append(out)
            print(f'IN: {a1}{a0}, out: {s1}{s0}')
    
    print()


# fig, ax = plt.subplots()
# im = ax.imshow(sums)
# ax.set_xticks(np.arange(len(x_axis)))
# ax.set_yticks(np.arange(len(y_axis)))
# ax.set_xticklabels(x_axis)
# ax.set_yticklabels(y_axis)
# plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
#          rotation_mode="anchor")
# for i in range(len(x_axis)):
#     for j in range(len(y_axis)):
#         text = ax.text(j, i, sums[i][j], ha="center", va="center", color="w")

# ax.set_title("CLA Adder Heat map")
# fig.tight_layout()
# plt.show()


# plt.figure()
# plt.plot(T,out)
# plt.show()