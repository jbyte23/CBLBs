from scipy.integrate import ode
import matplotlib.pyplot as plt
from alu_model import *
from models import *
from parameters import *
import numpy as np


def int_to_bool_list(num):
    return [int(bool(num & (1 << n))) for n in range(4)]


def list_to_int(xs):
    v = 0
    f = 1
    for x in xs:
        v += x*f
        f *= 2
    return v

rho_x = 0

params = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X

# simulation parameters
t_end = 1000
N = t_end


sums = []
x_axis = []
y_axis = []
for a in range(4):
    A = int_to_bool_list(a)
    A0 = A[0]
    A1 = A[1]

    sums.insert(0, [])
    x_axis.append(a)
    y_axis.insert(0, a)
    for b in range(4): 
        """
        ADDER

        A0, A1, B0, B1, \
        S0_I0, S0_I1, S0_J0, \
        S1_I0, S1_I1, S1_I2, S1_I3, S1_I4, S1_I5, S1_I6, S1_I7, S1_I8, \
        C_I0, C_I1, C_I2, C_I3, C_I4, C_I5, \
        S0_L_A, S0_L_B, S0_L_I0, S0_L_I1, S0_L_J0, \
        S1_L_A0_0, S1_L_A0_1, S1_L_A1_0, S1_L_A1_1, S1_L_B0_0, S1_L_B0_1, S1_L_B1_0, S1_L_B1_1, S1_L_I0, S1_L_I1, S1_L_I2, S1_L_I3, S1_L_I4, S1_L_I5, S1_L_I6, S1_L_I7, S1_L_I8, \
        C_L_A0_0, C_L_A0_1, C_L_B0_0, C_L_B0_1, C_L_I0, C_L_I1, C_L_I2, C_L_I3, C_L_I4, C_L_I5, \
        S0_N_A, S0_N_B, S0_N_I0, S0_N_I1, S0_N_J0, \
        S1_N_A0_0, S1_N_A0_1, S1_N_A1_0, S1_N_A1_1, S1_N_A1_2, S1_N_A1_3, S1_N_B0_0, S1_N_B0_1, S1_N_B1_0, S1_N_B1_1, S1_N_B1_2, S1_N_B1_3, S1_N_I0, S1_N_I1, S1_N_I2, S1_N_I3, S1_N_I4, S1_N_I5, S1_N_I6, S1_N_I7, S1_N_I8,\
        C_N_A0_0, C_N_A0_1, C_N_A1_0, C_N_A1_1, C_N_B0_0, C_N_B0_1, C_N_B1_0, C_N_B1_1, C_N_I0, C_N_I1, C_N_I2, C_N_I3, C_N_I4, C_N_I5, \
        S0, S1, C 
        """

        B = int_to_bool_list(b)
        B0 = B[0]
        B1 = B[1]
        

        Y0 = np.zeros(97)

        Y0[:4] = A0, A1, B0, B1
        # N_X nastavimo vse na 1 - st celic
        Y0[54:-3] = 1

        T = np.linspace(0, t_end, N)

        t1 = t_end
        dt = t_end/N
        T = np.arange(0,t1+dt,dt)
        Y = np.zeros([1+N,len(Y0)])
        Y[0,:] = Y0

        r = ode(adder_model_ODE).set_integrator('zvode', method='bdf')
        r.set_initial_value(Y0, T[0]).set_f_params(params)
        
        i = 1
        while r.successful() and r.t < t1:
            Y[i,:] = r.integrate(r.t+dt)
            i += 1


        S0 = 1 if Y[:,-3][-1] > 1 else 0
        S1 = 1 if Y[:,-2][-1] > 1 else 0
        C = 1 if Y[:,-1][-1] > 1 else 0

        s = C * 4 + S1 * 2 + S0

        sums[0].append(s)

        print(f' {A1}{A0}')
        print(f' {B1}{B0}')
        print('--------')
        print(f'{C}{S1}{S0}')
        print(f'{a}+{b}={s}')
        print()
    
    print()

print(x_axis, y_axis, sums)
fig, ax = plt.subplots()
plt.xlabel('A')
plt.ylabel('B')
im = ax.imshow(sums)
ax.set_xticks(np.arange(len(x_axis)))
ax.set_yticks(np.arange(len(y_axis)))
ax.set_xticklabels(x_axis)
ax.set_yticklabels(y_axis)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
for i in range(len(y_axis)):
    for j in range(len(x_axis)):
        text = ax.text(j, i, sums[i][j], ha="center", va="center", color="w")

ax.set_title("Adder Heat map")
fig.tight_layout()
_, labels = plt.yticks()
plt.show()
# plt.savefig('cblb/figs/heatmap_adder.pdf')

# plt.figure()
# plt.plot(T,out)
# plt.show()