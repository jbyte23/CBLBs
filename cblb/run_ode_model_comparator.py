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


comparisons = []
x_axis = []
y_axis = []
for a in range(4):
    A = int_to_bool_list(a)
    A0 = A[0]
    A1 = A[1]

    comparisons.insert(0, [])
    x_axis.append(a)
    y_axis.insert(0, a)
    for b in range(4): 
        """
        COMPARATOR

        A0, A1, B0, B1, \
        I0, I1, I0_0, I1_0,\
        L_A0, L_A1, L_B0, L_B1, L_I0, L_I1, L_I0_0, L_I1_0,\
        N_A0, N_A1, N_B0, N_B1, N_I0, N_I1, N_I0_0, N_I1_0,\
        S0, S1
        """

        B = int_to_bool_list(b)
        B0 = B[0]
        B1 = B[1]
        
        Y0 = np.zeros(26)

        Y0[:4] = A0, A1, B0, B1
        # N_X nastavimo vse na 1 - st celic
        Y0[16:-2] = 1

        T = np.linspace(0, t_end, N)

        t1 = t_end
        dt = t_end/N
        T = np.arange(0,t1+dt,dt)
        Y = np.zeros([1+N,len(Y0)])
        Y[0,:] = Y0

        r = ode(comparator_model_ODE).set_integrator('zvode', method='bdf')
        r.set_initial_value(Y0, T[0]).set_f_params(params)
        
        i = 1
        while r.successful() and r.t < t1:
            Y[i,:] = r.integrate(r.t+dt)
            i += 1

        S0 = 1 if Y[:,-2][-1] > 1 else 0
        S1 = 1 if Y[:,-1][-1] > 1 else 0

        s = S1 * 2 + S0

        comparisons[0].append(s)
        res = '='
        if s == 1:
            res = '<'
        elif s == 2:
            res = '>'

        print(f'{a} {res} {b}')
    
    print()

print(x_axis, y_axis, comparisons)
fig, ax = plt.subplots()
plt.xlabel('A')
plt.ylabel('B')
im = ax.imshow(comparisons)
ax.set_xticks(np.arange(len(x_axis)))
ax.set_yticks(np.arange(len(y_axis)))
ax.set_xticklabels(x_axis)
ax.set_yticklabels(y_axis)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
for i in range(len(y_axis)):
    for j in range(len(x_axis)):
        text = ax.text(j, i, comparisons[i][j], ha="center", va="center", color="w")

ax.set_title("Comparator Heat map")
fig.tight_layout()
_, labels = plt.yticks()
# plt.show()
plt.savefig('cblb/figs/heatmap_comparator.pdf')

# plt.figure()
# plt.plot(T,out)
# plt.show()