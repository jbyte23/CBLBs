import numpy as np


"""
YES and NOT cells
"""
def yes_cell(state, params):
    x, y, N_X, N_Y = state
    gamma_x, n_y, theta_x, delta_x, rho_x = params

    # presume that the molecules are degraded in the same strain as they are produced
    N_Y = N_X


    dx_dt = N_X * gamma_x * (y ** n_y)/(1 + (theta_x*y)**n_y ) - N_Y * (delta_x * x) - rho_x * x
    
    return dx_dt


def not_cell(state, params):
    L_X, x, y, N_X, N_Y = state
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x = params

    # presume that the molecules are degraded in the same strain as they are produced
    N_Y = N_X


    f = gamma_L_X * (y ** n_y)/(1 + (theta_L_X*y)**n_y )
    dL_X_dt = N_X * (f - delta_L * L_X)

    dx_dt = N_X * (eta_x * (1/(1+ (omega_x*L_X)**m_x))) - N_Y * (delta_x * x) - rho_x * x

    return dL_X_dt, dx_dt


# L_A ... intermediate
# a ... out
# b ... in
# N_A ... number of cells
def not_cell_wrapper(state, params):
    L_A, a, b, N_A = state
    
    state_A = L_A, a, b, N_A, N_A
    params_A = params

    return not_cell(state_A, params_A)

# a ... out
# b ... in
# N_A ... number of cells
def yes_cell_wrapper(state, params):
    a, b, N_A = state

    state_A = a, b, N_A, N_A
    params_A = params

    return yes_cell(state_A, params_A)


"""
OR cells 2 input
"""

def yes_yes_or2(state, params):
    # YES x OR YES y = OR
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X = params
    params_yes = gamma_x, n_y, theta_x, delta_x, rho_x

    A0, B0, N_A0, N_B0, OR_out = state

    state_A = OR_out, A0, N_A0
    state_B = OR_out, B0, N_B0

    dOR_out = 0
    dOR_out += yes_cell_wrapper(state_A, params_yes)
    dN_A0 = 0
    dOR_out += yes_cell_wrapper(state_B, params_yes)
    dN_B0 = 0
            
    return dN_A0, dN_B0, dOR_out


def yes_not_or2(state, params):
    # YES x OR YES y = OR
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X = params
    params_yes = gamma_x, n_y, theta_x, delta_x, rho_x
    params_not = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x

    A0, B0, L_B0, N_A0, N_B0, OR_out = state

    state_A = OR_out, A0, N_A0
    state_B = L_B0, OR_out, B0, N_B0

    dOR_out = 0
    dOR_out += yes_cell_wrapper(state_A, params_yes)
    dN_A0 = 0
    dL_B0, dd = not_cell_wrapper(state_B, params_not)
    dOR_out += dd
    dN_B0 = 0
            
    return dL_B0, dN_A0, dN_B0, dOR_out


def not_not_or2(state, params):
    # NAND
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X = params
    params_not = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x

    A, B, L_A, L_B, N_A, N_B, OR_out = state

    state_A = L_A, OR_out, A, N_A
    state_B = L_B, OR_out, B, N_B

    dOR_out = 0
    dL_A, dd = not_cell_wrapper(state_A, params_not)
    dOR_out += dd
    dN_A = 0
    dL_B, dd = not_cell_wrapper(state_B, params_not)
    dOR_out += dd
    dN_B = 0

    return dL_A, dL_B, dN_A, dN_B, dOR_out



"""
ALU components cells
"""
def adder_cell(state, params):
    return None


def shifter_cell(state, params):

    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X = params
    params_not = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x

    # Vhodi
    A0, A1, M0, M1, I0, I1, \
    L_A0, L_A1, L_M0, L_M1, L_I0, L_I1, \
    N_A0, N_A1, N_M0, N_M1, N_I0, N_I1, \
    S0, S1 = state


    """
    I0: NOT M0 OR NOT A1
    """
    state_I0 = A1, M0, L_A1, L_M0, N_A1, N_M0, I0
    I0_out = not_not_or2(state_I0, params)
    dL_A1, dL_M0, dN_A1, dN_M0, dI0 = I0_out

    """
    I1: NOT M1 OR NOT A0
    """
    state_I1 = A0, M1, L_A0, L_M1, N_A0, N_M1, I1
    I1_out = not_not_or2(state_I1, params)
    dL_A0, dL_M1, dN_A0, dN_M1, dI1 = I1_out
    

    """
    S0: NOT I0
    """
    state_S0 = L_I0, S0, I0, N_I0
    dL_I0, S0_out = not_cell_wrapper(state_S0, params_not)
    dN_I0 = 0

    """
    S1: NOT I1
    """
    state_S1 = L_I1, S1, I1, N_I1
    dL_I1, S1_out = not_cell_wrapper(state_S1, params_not)
    dN_I1 = 0

    dA0, dA1, dM0, dM1 = 0, 0, 0, 0

    return dA0, dA1, dM0, dM1, dI0, dI1, \
            dL_A0, dL_A1, dL_M0, dL_M1, dL_I0, dL_I1, \
            dN_A0, dN_A1, dN_M0, dN_M1, dN_I0, dN_I1, \
            S0_out, S1_out



def comparator_cell(state, params):
    return None




"""
ALU component models
"""

def adder_model(state, T, params):
    return None


def shifter_model(state, T, params):
    return np.array(shifter_cell(state, params))


def comparator_model(state, T, params):
    return None


"""
ALU model
"""

def alu_model(state, T, params):
    return None



"""
testing models
"""

def not_not_or2_model(state, T, params):

    # Vhodi
    # A0, A1, B0, B1 = state[:4]
    A0, B0, L_A0, L_B0, N_A0, N_B0, OR_out = state

    state_OR = A0, B0, L_A0, L_B0, N_A0, N_B0, OR_out
    dL_A0, dL_B0, dN_A0, dN_B0, dOR_out = not_not_or2(state_OR, params)

    dA0, dB0 = 0, 0

    return np.array([dA0, dB0,
                    dL_A0, dL_B0,
                    dN_A0, dN_B0,
                    dOR_out
    ])

def yes_yes_or2_model(state, T, params):

    # Vhodi
    A0, B0, N_A0, N_B0, OR_out = state

    state_OR0 = A0, B0, N_A0, N_B0, OR_out
    dN_A0, dN_B0, dOR_out = yes_yes_or2(state_OR0, params)

    dA0, dB0 = 0, 0

    return np.array([dA0, dB0,
                    dN_A0, dN_B0,
                    dOR_out
    ])

def yes_not_or2_model(state, T, params):

    # Vhodi
    # A0, A1, B0, B1 = state[:4]
    A0, B0, L_B0, N_A0, N_B0, OR_out = state

    state_OR = A0, B0, L_B0, N_A0, N_B0, OR_out
    dL_B0, dN_A0, dN_B0, dOR_out = yes_not_or2(state_OR, params)

    dA0, dB0 = 0, 0

    return np.array([dA0, dB0,
                    dL_B0,
                    dN_A0, dN_B0,
                    dOR_out
    ])

def two_bit_not_not_model(state, T, params):

    # Vhodi
    A0, A1, B0, B1, \
    L_A0, L_A1, L_B0, L_B1, \
    N_A0, N_A1, N_B0, N_B1, \
    S0, S1 = state

    state_S0 = A0, B0, L_A0, L_B0, N_A0, N_B0, S0
    state_S1 = A1, B1, L_A1, L_B1, N_A1, N_B1, S1

    S0_out = not_not_or2(state_S0, params)
    dL_A0, dL_B0, dN_A0, dN_B0, dS0 = S0_out
    S1_out = not_not_or2(state_S1, params)
    dL_A1, dL_B1, dN_A1, dN_B1, dS1 = S1_out

    dA0, dA1, dB0, dB1 = 0, 0, 0, 0

    return np.array([dA0, dA1, dB0, dB1, 
    dL_A0, dL_A1, dL_B0, dL_B1, 
    dN_A0, dN_A1, dN_B0, dN_B1, 
    dS0, dS1])

"""
Wrappers for scipy.integrate.ode
"""

def alu_model_ODE(T, state, params):
    return alu_model(state, T, params)

def adder_model_ODE(T, state, params):
    return adder_model(state, T, params)

def shifter_model_ODE(T, state, params):
    return shifter_model(state, T, params)

def comparator_model_ODE(T, state, params):
    return comparator_model(state, T, params)

def test_model_ODE(T, state, params):
    return yes_not_or2_model(state, T, params)