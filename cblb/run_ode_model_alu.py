from typing import Counter
from scipy.integrate import ode
import matplotlib.pyplot as plt
from alu_model import *
from models import *
from parameters import *
import numpy as np


def int_to_bool_list(num):
    return [int(bool(num & (1 << n))) for n in range(2)]


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
t_end = 400
N = t_end


sums = []
x_axis = []
y_axis = []
M0_out = []
M1_out = []
A0_out = []
A1_out = []
B0_out = []
B1_out = []
S0_out = []
S1_out = []
C_out = []

for m in range(4):
    
    if m == 0: print('#####\nSUM\n#####')
    elif m == 1: print('#####\nSHIFT RIGHT\n#####')
    elif m == 2: print('#####\nSHIFT LEFT\n#####')
    elif m == 3: print('#####\nCOMPARE\n#####')

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

            Y0 = np.zeros(203)

            M0, M1 = int_to_bool_list(m)
            Y0[:6] = A0, A1, B0, B1, M0, M1
            # N_X nastavimo vse na 1 - st celic
            Y0[56:96] = 1
            Y0[111:119] = 1
            Y0[129:135] = 1
            Y0[153:169] = 1
            Y0[186:202] = 1

            T = np.linspace(0, t_end, N)

            t1 = t_end
            dt = t_end/N
            T = np.arange(0, t1+dt, dt)
            Y = np.zeros([1+N, len(Y0)])
            Y[0, :] = Y0

            r = ode(alu_model_ODE).set_integrator('zvode', method='bdf')
            r.set_initial_value(Y0, T[0]).set_f_params(params)

            i = 1
            while r.successful() and r.t < t1:
                Y[i, :] = r.integrate(r.t+dt)
                i += 1

            S0 = 1 if Y[:, 169][-1] > 1 else 0
            S1 = 1 if Y[:, -1][-1] > 1 else 0
            C = 1 if Y[:, 98][-1] > 1 else 0


            if A0 == 1 and A1 == 1 and B0 == 0 and B1 == 1:
                A0_out.append(Y[:, 0])
                A1_out.append(Y[:, 1])
                B0_out.append(Y[:, 2])
                B1_out.append(Y[:, 3])
                S0_out.append(Y[:, 169])
                S1_out.append(Y[:, -1])
                C_out.append(Y[:, 98])

            s = C * 4 + S1 * 2 + S0

            sums[0].append(s)

            print(f' {A1}{A0}')
            print(f' {B1}{B0}')
            print(f'{C}{S1}{S0}')
            # print(f'{a}+{b}={s}')
            print()

        print()

    M0_out.append(Y[:, 4])
    M1_out.append( Y[:, 5])

# print(x_axis, y_axis, sums)
# fig, ax = plt.subplots()
# plt.xlabel('A')
# plt.ylabel('B')
# im = ax.imshow(sums)
# ax.set_xticks(np.arange(len(x_axis)))
# ax.set_yticks(np.arange(len(y_axis)))
# ax.set_xticklabels(x_axis)
# ax.set_yticklabels(y_axis)
# plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
#          rotation_mode="anchor")
# for i in range(len(y_axis)):
#     for j in range(len(x_axis)):
#         text = ax.text(j, i, sums[i][j], ha="center", va="center", color="w")

# ax.set_title("Adder Heat map")
# fig.tight_layout()
# _, labels = plt.yticks()
# plt.show()




# A0, A1= Y[:, 0], Y[:, 1]
# B0, B1 = Y[:, 2], Y[:, 3]

# M0, M1 = Y[:, 4], Y[:, 5]

# S0 = Y[:, 169]
# S1 = Y[:, -1]
C = Y[:, 98]

print(len(T))
# plot
plt.figure(1)
plt.title('Full adder')
ax2 = plt.subplot(711)
ax2.plot(T, M0_out[0], color="#ff1a1a", alpha=0.7)
ax2.plot(T, M1_out[0], color="#3366ff", alpha=0.7)
ax2.legend(["$M_0$", "$M_1$"])

ax1 = plt.subplot(712)
ax1.plot(T, A0_out[0], color="#ff1a1a", alpha=0.7)
ax1.plot(T, A1_out[0], color="#3366ff", alpha=0.7)
ax1.legend(["$A_0$", "$A_1$"])

ax5 = plt.subplot(713)
ax5.plot(T, B0_out[0], color="#ff1a1a", alpha=0.7)
ax5.plot(T, B1_out[0], color="#3366ff", alpha=0.7)
ax5.legend(["$B_0$", "$B_1$"])

ax6 = plt.subplot(714)
ax6.plot(T, C_out[0], color="#663300", alpha=0.7)
ax6.legend(["$C_{OUT}$"])

ax9 = plt.subplot(715)
ax9.plot(T, S1_out[0], color="#663300", alpha=0.7)
ax9.legend(["$S_1$"])

ax10 = plt.subplot(716)
ax10.plot(T, S0_out[0], color="#663300", alpha=0.7)
ax10.legend(["$S_0$"])

text_params = {'ha': 'center', 'va': 'center', 'family': 'sans-serif', 'fontweight': 'bold'}
# plt.text(200, 30, '$A_0=' + str(int(A0[0])) + '$  $A_1=' + str(int(A1[0])), color='r', size=20, **text_params)
# plt.text(200, 25, '$B_0=' + str(int(B0[0])) + '$  $B_1=' + str(int(B1[0])), color='r', size=20, **text_params)

ticks = np.arange(0, 12, 1) * t_end
tick_labels = list(map(str, ticks))
for i in range(len(tick_labels)):
    if i % 2 == 1:
        tick_labels[i] = ""

plt.gcf().set_size_inches(8, 10)
plt.savefig("cblb/figs/graphs_adder.pdf", bbox_inches='tight')

# plt.show()



# plot
plt.figure(2)
plt.title('Shift right')
ax2 = plt.subplot(711)
ax2.plot(T, M0_out[1], color="#ff1a1a", alpha=0.7)
ax2.plot(T, M1_out[1], color="#3366ff", alpha=0.7)
ax2.legend(["$M_0$", "$M_1$"])

ax1 = plt.subplot(712)
ax1.plot(T, A0_out[1], color="#ff1a1a", alpha=0.7)
ax1.plot(T, A1_out[1], color="#3366ff", alpha=0.7)
ax1.legend(["$A_0$", "$A_1$"])

ax9 = plt.subplot(713)
ax9.plot(T, S1_out[1], color="#663300", alpha=0.7)
ax9.legend(["$S_1$"])

ax10 = plt.subplot(714)
ax10.plot(T, S0_out[1], color="#663300", alpha=0.7)
ax10.legend(["$S_0$"])

text_params = {'ha': 'center', 'va': 'center', 'family': 'sans-serif', 'fontweight': 'bold'}
# plt.text(200, 30, '$A_0=' + str(int(A0[0])) + '$  $A_1=' + str(int(A1[0])), color='r', size=20, **text_params)
# plt.text(200, 25, '$B_0=' + str(int(B0[0])) + '$  $B_1=' + str(int(B1[0])), color='r', size=20, **text_params)

ticks = np.arange(0, 12, 1) * t_end
tick_labels = list(map(str, ticks))
for i in range(len(tick_labels)):
    if i % 2 == 1:
        tick_labels[i] = ""

plt.gcf().set_size_inches(8, 10)
plt.savefig("cblb/figs/graph_shift_left.pdf", bbox_inches='tight')

# plt.show()



# plot
plt.figure(3)
plt.title('Shift left')
ax2 = plt.subplot(711)
ax2.plot(T, M0_out[2], color="#ff1a1a", alpha=0.7)
ax2.plot(T, M1_out[2], color="#3366ff", alpha=0.7)
ax2.legend(["$M_0$", "$M_1$"])

ax1 = plt.subplot(712)
ax1.plot(T, A0_out[2], color="#ff1a1a", alpha=0.7)
ax1.plot(T, A1_out[2], color="#3366ff", alpha=0.7)
ax1.legend(["$A_0$", "$A_1$"])

ax9 = plt.subplot(713)
ax9.plot(T, S1_out[2], color="#663300", alpha=0.7)
ax9.legend(["$S_1$"])

ax10 = plt.subplot(714)
ax10.plot(T, S0_out[2], color="#663300", alpha=0.7)
ax10.legend(["$S_0$"])

text_params = {'ha': 'center', 'va': 'center', 'family': 'sans-serif', 'fontweight': 'bold'}
# plt.text(200, 30, '$A_0=' + str(int(A0[0])) + '$  $A_1=' + str(int(A1[0])), color='r', size=20, **text_params)
# plt.text(200, 25, '$B_0=' + str(int(B0[0])) + '$  $B_1=' + str(int(B1[0])), color='r', size=20, **text_params)

ticks = np.arange(0, 12, 1) * t_end
tick_labels = list(map(str, ticks))
for i in range(len(tick_labels)):
    if i % 2 == 1:
        tick_labels[i] = ""

plt.gcf().set_size_inches(8, 10)
plt.savefig("cblb/figs/graph_shift_right.pdf", bbox_inches='tight')

# plt.show()



# plot
plt.figure(4)
plt.title('Comparator')
ax2 = plt.subplot(711)
ax2.plot(T, M0_out[3], color="#ff1a1a", alpha=0.7)
ax2.plot(T, M1_out[3], color="#3366ff", alpha=0.7)
ax2.legend(["$M_0$", "$M_1$"])

ax1 = plt.subplot(712)
ax1.plot(T, A0_out[3], color="#ff1a1a", alpha=0.7)
ax1.plot(T, A1_out[3], color="#3366ff", alpha=0.7)
ax1.legend(["$A_0$", "$A_1$"])

ax5 = plt.subplot(713)
ax5.plot(T, B0_out[3], color="#ff1a1a", alpha=0.7)
ax5.plot(T, B1_out[3], color="#3366ff", alpha=0.7)
ax5.legend(["$B_0$", "$B_1$"])

ax9 = plt.subplot(714)
ax9.plot(T, S1_out[3], color="#663300", alpha=0.7)
ax9.legend(["$S_1$"])

ax10 = plt.subplot(715)
ax10.plot(T, S0_out[3], color="#663300", alpha=0.7)
ax10.legend(["$S_0$"])

text_params = {'ha': 'center', 'va': 'center', 'family': 'sans-serif', 'fontweight': 'bold'}
# plt.text(200, 30, '$A_0=' + str(int(A0[0])) + '$  $A_1=' + str(int(A1[0])), color='r', size=20, **text_params)
# plt.text(200, 25, '$B_0=' + str(int(B0[0])) + '$  $B_1=' + str(int(B1[0])), color='r', size=20, **text_params)

ticks = np.arange(0, 12, 1) * t_end
tick_labels = list(map(str, ticks))
for i in range(len(tick_labels)):
    if i % 2 == 1:
        tick_labels[i] = ""

plt.gcf().set_size_inches(8, 10)
plt.savefig("cblb/figs/graph_comparator.pdf", bbox_inches='tight')

# plt.show()

