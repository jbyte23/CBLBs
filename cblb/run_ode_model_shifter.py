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


shifts = []
x_axis = list(range(4))
y_axis = [2,1]
for m in [1,2]: # 1 (01): right, 2 (10): left
    
    M = int_to_bool_list(m)
    M0 = M[0]
    M1 = M[1]
    print(f'MODE: {M1}{M0}')

    shifts.append([])
    # x_axis.append(m)
    # y_axis.insert(0, m)
    for a in range(4): 
        """
        SHIFTER

        A0, A1, M0, M1, I0, I1, \
        L_A0, L_A1, L_M0, L_M1, L_I0, L_I1, \
        N_A0, N_A1, N_M0, N_M1, N_I0, N_I1, \
        S0, S1
        """

        A = int_to_bool_list(a)
        A0 = A[0]
        A1 = A[1]
        
        Y0 = np.zeros(20)

        Y0[:4] = A0, A1, M0, M1
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

        S0 = 1 if Y[:,-2][-1] > 2 else 0
        S1 = 1 if Y[:,-1][-1] > 2 else 0

        s = S1 * 2 + S0

        shifts[-1].append(s)
        print(f'IN: {A1}{A0}, out: {S1}{S0}')
    
    print()

print(x_axis, y_axis, shifts)
fig, ax = plt.subplots()
im = ax.imshow(shifts)
ax.set_xticks(np.arange(len(x_axis)))
ax.set_yticks(np.arange(len(y_axis)))
ax.set_xticklabels(x_axis)
ax.set_yticklabels(y_axis)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
for i in range(len(y_axis)):
    for j in range(len(x_axis)):
        text = ax.text(j, i, shifts[i][j], ha="center", va="center", color="w")

ax.set_title("Shifter Heat map")
fig.tight_layout()
plt.show()


# plt.figure()
# plt.plot(T,out)
# plt.show()