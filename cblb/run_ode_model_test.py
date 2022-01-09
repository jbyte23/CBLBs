from scipy.integrate import ode
import matplotlib.pyplot as plt
from alu_model import *
from models import *
from parameters import *
import numpy as np

rho_x = 0

params = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X

# Vhodne logicne vrednosti
A = [0,1]
B = [0,1]

"""
yes_not_or_2
"""
# # Y = A0, B0, L_B0, N_A0, N_B0, OR_out
# Y0 = np.zeros(6)

# Y0[:2] = A, B
# # N_X nastavimo vse na 1 - st celic
# Y0[2:] = 1


"""
yes_yes_or_2
"""
# Y = A0, B0, N_A0, N_B0, OR_out
# Y0 = np.zeros(5)

# Y0[:2] = A, B
# # N_X nastavimo vse na 1 - st celic
# Y0[2:] = 1




# simulation parameters
t_end = 1000
N = t_end


sums = []
x_axis = []
y_axis = []

for a in A:
    sums.insert(0, [])
    x_axis.append(a)    
    y_axis.insert(0, a)

    for b in B:
        """
        not_not_or_2
        """
        # # Y = A, B, L_A, L_B, N_A, N_B, OR_out
        # Y0 = np.zeros(7)

        # Y0[:2] = a, b
        # # N_X nastavimo vse na 1 - st celic
        # Y0[4:6] = 1

        """
        yes_not_or_2
        """
        # # Y = A0, B0, L_B0, N_A0, N_B0, OR_out
        # Y0 = np.zeros(6)

        # Y0[:2] = a, b
        # # N_X nastavimo vse na 1 - st celic
        # Y0[3:5] = 1


        """
        yes_yes_or_2
        """
        # # Y = A0, B0, N_A0, N_B0, OR_out
        # Y0 = np.zeros(5)

        # Y0[:2] = a, b
        # # N_X nastavimo vse na 1 - st celic
        # Y0[2:] = 1

        """"
        not_not_not_or3

        A, B, C, \
        L_A, L_B, L_C, \
        N_A, N_B, N_C, \
        OR_out
        """
        # Y0 = np.zeros(10)
        # c = 1
        # Y0[:3] = a, b, c
        # # N_X nastavimo vse na 1 - st celic
        # Y0[6:10] = 1


        """
        not_not_not_yes_or4

        A, B, C, D, \
        L_A, L_B, L_C, \
        N_A, N_B, N_C, N_D, \
        OR_out
        """
        # Y0 = np.zeros(12)
        # c = 1
        # d = 0
        # Y0[:4] = a, b, c, d
        # # N_X nastavimo vse na 1 - st celic
        # Y0[7:11] = 1

        """
        not_not_not_not_or4

        A, B, C, D, \
        L_A, L_B, L_C, L_D, \
        N_A, N_B, N_C, N_D, \
        OR_out
        """
        Y0 = np.zeros(13)
        c = 1
        d = 1
        Y0[:4] = a, b, c, d
        # N_X nastavimo vse na 1 - st celic
        Y0[5:12] = 1

        """
        two_bit_not_not_or_2

        A0, A1, B0, B1, \
        L_A0, L_A1, L_B0, L_B1, \
        N_A0, N_A1, N_B0, N_B1, \
        S0, S1
        """
        # Y0 = np.zeros(14)

        # a0 = a
        # b0 = 1  # se negira
        # a1 = 0  # se negira
        # b1 = b

        # Y0[:4] = a0, a1, b0, b1
        # Y0[8:12] = 1




        T = np.linspace(0, t_end, N)

        t1 = t_end
        dt = t_end/N
        T = np.arange(0,t1+dt,dt)
        Y = np.zeros([1+N,len(Y0)])
        Y[0,:] = Y0

        r = ode(test_model_ODE).set_integrator('zvode', method='bdf')
        r.set_initial_value(Y0, T[0]).set_f_params(params)

        i = 1
        while r.successful() and r.t < t1:
            Y[i,:] = r.integrate(r.t+dt)
            i += 1

        # 1 bit
        out = Y[:,-1][-1]
        sums[0].append(out > 1)
        # or2
        # print(f'in: {a}{b}, out: {out > 1} ({out})')
        # or3
        # print(f'in: {a}{b}{c}, out: {out > 1} ({out})')
        # or4
        print(f'in: {a}{b}{c}{d}, out: {out > 1} ({out})')

        # 2 bit
        # S0 = Y[:,-2][-1] > 1
        # S1 = Y[:,-1][-1] > 1
        # print(f'a0b0: {int(Y0[0])}{int(Y0[2])}, out: {S0}\na1b1: {int(Y0[1])}{int(Y0[3])}, out: {S1}')
        # print()


fig, ax = plt.subplots()
im = ax.imshow(sums)
ax.set_xticks(np.arange(len(x_axis)))
ax.set_yticks(np.arange(len(y_axis)))
ax.set_xticklabels(x_axis)
ax.set_yticklabels(y_axis)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
for i in range(len(x_axis)):
    for j in range(len(y_axis)):
        text = ax.text(j, i, sums[i][j], ha="center", va="center", color="w")

ax.set_title("CLA Adder Heat map")
fig.tight_layout()
plt.show()


# plt.figure()
# plt.plot(T,out)
# plt.show()