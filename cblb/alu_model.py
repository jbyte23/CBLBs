import numpy as np

def not_cell(state, params):
    L_X, x, y = state
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, N_X, N_Y = params

    # presume that the molecules are degraded in the same strain as they are produced
    N_Y = N_X

    f = gamma_L_X * (y ** n_y)/(1 + (theta_L_X*y)**n_y )
    dL_X_dt = N_X * (f - delta_L * L_X)

    dx_dt = N_X * (eta_x * (1/(1+ (omega_x*L_X)**m_x))) - N_Y * (delta_x * x) - rho_x * x

    return dL_X_dt, dx_dt

def yes_cell(state, params):
    x, y = state
    gamma_x, n_y, theta_x, delta_x, rho_x, N_X, N_Y = params

    # presume that the molecules are degraded in the same strain as they are produced
    N_Y = N_X

    dx_dt = N_X * gamma_x * (y ** n_y)/(1 + (theta_x*y)**n_y ) - N_Y * (delta_x * x) - rho_x * x
    
    return dx_dt

# L_A ... intermediate
# a ... out
# b ... in
def not_cell_wrapper(state, params):
    L_A, a, b = state
    
    state_A = L_A, a, b
    
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, N_X = params
    params_A = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, N_X, N_X
    

    return not_cell(state_A, params_A)

# a ... out
# b ... in
def yes_cell_wrapper(state, params):
    a, b = state

    state_A = a, b
    gamma_x, n_y, theta_x, delta_x, rho_x, N_X = params
    params_A = gamma_x, n_y, theta_x, delta_x, rho_x, N_X, N_X
    

    return yes_cell(state_A, params_A)

def not_model(state, T, params):
    L_A, a, b = state

    delta_L, gamma_L_A, n_b, theta_L_A, eta_a, omega_a, m_a, delta_a, delta_b, rho_a, rho_b, N_A = params

    state_not = L_A, a, b
    params_not = delta_L, gamma_L_A, n_b, theta_L_A, eta_a, omega_a, m_a, delta_a, rho_a, N_A
    dL_A_dt, da_dt = not_cell_wrapper(state_not, params_not)
    
    db_dt = 0#- N_A * delta_b * b - rho_b * b

    return np.array([dL_A_dt, da_dt, db_dt])

def yes_model(state, T, params):
    a, b = state
    
    gamma_a, n_b, theta_a, delta_a, delta_b, rho_a, rho_b, N_A = params

    state_yes = a, b
    params_yes = gamma_a, n_b, theta_a, delta_a, rho_a, N_A
    da_dt = yes_cell_wrapper(state_yes, params_yes)
    
    db_dt = 0 #- N_A * delta_b * b - rho_b * b

    return np.array([da_dt, db_dt])



def yes_yes_or2(state, params):
    # YES x OR YES y = OR

    OR_out, A0, N_A0, B0, N_B0 = state

    state_A = OR_out, A0, N_A0
    state_B = OR_out, B0, N_B0

    dOR_out = 0
    dOR_out += yes_cell_wrapper(state_A, params)
    dN_A0 = 0
    dOR_out += yes_cell_wrapper(state_B, params)
    dN_B0 = 0
            

    return None

def or4(state, params):
    return None

def adder_model(state, T, params):
    return None

def shifter_model(state, T, params):
    return None

def comparator_model(state, T, params):
    return None



def ale_model(state, T, params):
    return None

def test_model(state, T, params):

    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X = params
    params_yes = gamma_x, n_y, theta_x, delta_x, rho_x
    # params_not = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x

    # Vhodi
    # A0, A1, B0, B1 = state[:4]
    A0, B0 = state[:2]
    # Parameter za skupek celic (npr. I0)
    OR_out = state[2]
    N_A0, N_B0 = state[3:5]
    # out = state[5]

    state = OR_out, A0, N_A0, B0, N_B0
    yes_yes_or = yes_yes_or2(state, params)


"""
Wrappers for scipy.integrate.ode
"""

def ale_model_ODE(T, state, params):
    return ale_model(state, T, params)

def adder_model_ODE(T, state, params):
    return adder_model(state, T, params)

def shifter_model_ODE(T, state, params):
    return shifter_model(state, T, params)

def comparator_model_ODE(T, state, params):
    return comparator_model(state, T, params)