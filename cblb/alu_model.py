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
AUXILIARY CELLS
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


def not_not_not_or3(state, params):
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X = params
    params_yes = gamma_x, n_y, theta_x, delta_x, rho_x
    params_not = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x

    A, B, C, \
    L_A, L_B, L_C, \
    N_A, N_B, N_C, \
    OR_out = state

    dOR_out = 0

    state_A = L_A, OR_out, A, N_A
    dL_A, dd = not_cell_wrapper(state_A, params_not)
    dOR_out += dd
    dN_A = 0

    state_B = L_B, OR_out, B, N_B
    dL_B, dd = not_cell_wrapper(state_B, params_not)
    dOR_out += dd
    dN_B = 0

    state_C = L_C, OR_out, C, N_C
    dL_C, dd = not_cell_wrapper(state_C, params_not)
    dOR_out += dd
    dN_C = 0

    return dL_A, dL_B, dL_C, dN_A, dN_B, dN_C, dOR_out


def not_not_not_yes_or4(state, params):

    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X = params
    params_yes = gamma_x, n_y, theta_x, delta_x, rho_x
    params_not = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x

    A, B, C, D, \
    L_A, L_B, L_C, \
    N_A, N_B, N_C, N_D, \
    OR_out = state

    dOR_out = 0

    state_A = L_A, OR_out, A, N_A
    dL_A, dd = not_cell_wrapper(state_A, params_not)
    dOR_out += dd
    dN_A = 0

    state_B = L_B, OR_out, B, N_B
    dL_B, dd = not_cell_wrapper(state_B, params_not)
    dOR_out += dd
    dN_B = 0

    state_C = L_C, OR_out, C, N_C
    dL_C, dd = not_cell_wrapper(state_C, params_not)
    dOR_out += dd
    dN_C = 0

    state_D = OR_out, D, N_D
    dOR_out += yes_cell_wrapper(state_D, params_yes)
    dN_D = 0

    return dL_A, dL_B, dL_C, dN_A, dN_B, dN_C, dN_D, dOR_out


def not_not_not_not_or4(state, params):

    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X = params
    params_not = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x

    A, B, C, D, \
    L_A, L_B, L_C, L_D, \
    N_A, N_B, N_C, N_D, \
    OR_out = state

    dOR_out = 0

    state_A = L_A, OR_out, A, N_A
    dL_A, dd = not_cell_wrapper(state_A, params_not)
    dOR_out += dd
    dN_A = 0

    state_B = L_B, OR_out, B, N_B
    dL_B, dd = not_cell_wrapper(state_B, params_not)
    dOR_out += dd
    dN_B = 0

    state_C = L_C, OR_out, C, N_C
    dL_C, dd = not_cell_wrapper(state_C, params_not)
    dOR_out += dd
    dN_C = 0

    state_D = L_D, OR_out, D, N_D
    dL_D, dd = not_cell_wrapper(state_D, params_not)
    dOR_out += dd
    dN_D = 0

    return dL_A, dL_B, dL_C, dL_D, dN_A, dN_B, dN_C, dN_D, dOR_out


def sum0(state, params):
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X = params
    params_not = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x

    A, B, \
    I0, I1, J0, \
    L_A, L_B, L_I0, L_I1, L_J0, \
    N_A, N_B, N_I0, N_I1, N_J0, \
    S0 = state

    """
    I0: NOT A0 OR NOT B0
    """
    state_I0 = A, B, L_A, L_B, N_A, N_B, I0
    dL_A, dL_B, dN_A, dN_B, dI0 = not_not_or2(state_I0, params)

    """
    I1: A0 OR B0
    """
    state_I1 = A, B, N_A, N_B, I1
    dN_A, dN_B, dI1 = yes_yes_or2(state_I1, params)

    """
    J0: NOT I0 OR NOT I1
    """
    state_J0 = I0, I1, L_I0, L_I1, N_I0, N_I1, J0
    dL_I0, dL_I1, dN_I0, dN_I1, dJ0 = not_not_or2(state_J0, params)


    """
    S0: NOT J0
    """
    state_S0 = L_J0, S0, J0, N_J0
    dL_J0, dS0 = not_cell_wrapper(state_S0, params_not)
    dN_J0 = 0


    dA, dB = 0, 0,

    return dA, dB, \
            dI0, dI1, dJ0, \
            dL_A, dL_B, dL_I0, dL_I1, dL_J0, \
            dN_A, dN_B, dN_I0, dN_I1, dN_J0, \
            dS0


def carry_out(state, params):
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X = params
    params_not = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x

    A0, A1, B0, B1, \
    I0, I1, I2, I3, I4, I5, \
    L_A0_0, L_A0_1, L_B0_0, L_B0_1, L_I0, L_I1, L_I2, L_I3, L_I4, L_I5, \
    N_A0_0, N_A0_1, N_A1_0, N_A1_1, N_B0_0, N_B0_1, N_B1_0, N_B1_1, N_I0, N_I1, N_I2, N_I3, N_I4, N_I5, \
    C = state

    """
    I0: NOT A0 NOT B0
    """
    state_I0 = A0, B0, L_A0_0, L_B0_0, N_A0_0, N_B0_0, I0
    dL_A0_0, dL_B0_0, dN_A0_0, dN_B0_0, dI0 = not_not_or2(state_I0, params)

    """
    I1: NOT A0 NOT B0
    """
    state_I1 = A0, B0, L_A0_1, L_B0_1, N_A0_1, N_B0_1, I1
    dL_A0_1, dL_B0_1, dN_A0_1, dN_B0_1, dI1 = not_not_or2(state_I1, params)

    """
    I2: A1 OR B1
    """
    state_I2 = A1, B1, N_A1_0, N_B1_0, I2
    dN_A1_0, dN_B1_0, dI2 = yes_yes_or2(state_I2, params)

    """
    I3: A1 OR NOT I0
    """
    state_I3 = A1, I0, L_I0, N_A1_1, N_I0, I3
    dL_I0, dN_A1_1, dN_I0, dI3 = yes_not_or2(state_I3, params)

    """
    I4: B1 OR NOT I1
    """
    state_I4 = B1, I1, L_I1, N_B1_1, N_I1, I4
    dL_I1, dN_B1_1, dN_I1, dI4 = yes_not_or2(state_I4, params)

    """
    I5: NOT I2 OR NOT I3 OR NOT I4
    """
    state_I5 = I2, I3, I4, L_I2, L_I3, L_I4, N_I2, N_I3, N_I4, I5
    dL_I2, dL_I3, dL_I4, dN_I2, dN_I3, dN_I4, dI5 = not_not_not_or3(state_I5, params)


    """
    C: NOT I5
    """
    state_C = L_I5, C, I5, N_I5
    dL_I5, dC = not_cell_wrapper(state_C, params_not)
    dN_I5 = 0

    dA0, dA1, dB0, dB1 = 0, 0, 0, 0

    return dA0, dA1, dB0, dB1, \
            dI0, dI1, dI2, dI3, dI4, dI5, \
            dL_A0_0, dL_A0_1, dL_B0_0, dL_B0_1, dL_I0, dL_I1, dL_I2, dL_I3, dL_I4, dL_I5, \
            dN_A0_0, dN_A0_1, dN_A1_0, dN_A1_1, dN_B0_0, dN_B0_1, dN_B1_0, dN_B1_1, dN_I0, dN_I1, dN_I2, dN_I3, dN_I4, dN_I5, \
            dC


def sum1(state, params):
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X = params
    params_not = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x

    A0, A1, B0, B1, \
    I0, I1, I2, I3, I4, I5, I6, I7, I8, \
    L_A0_0, L_A0_1, L_A1_0, L_A1_1, L_B0_0, L_B0_1, L_B1_0, L_B1_1, L_I0, L_I1, L_I2, L_I3, L_I4, L_I5, L_I6, L_I7, L_I8, \
    N_A0_0, N_A0_1, N_A1_0, N_A1_1, N_A1_2, N_A1_3, N_B0_0, N_B0_1, N_B1_0, N_B1_1, N_B1_2, N_B1_3, N_I0, N_I1, N_I2, N_I3, N_I4, N_I5, N_I6, N_I7, N_I8,\
    S1 = state

    """
    I0: NOT A1 OR NOT B1
    """
    state_I0 = A1, B1, L_A1_0, L_B1_0, N_A1_0, N_B1_0, I0
    dL_A1_0, dL_B1_0, dN_A1_0, dN_B1_0, dI0 = not_not_or2(state_I0, params)

    """
    I1: A1 OR B1
    """
    state_I1 = A1, B1, N_A1_1, N_B1_1, I1
    dN_A1_1, dN_B1_1, dI1 = yes_yes_or2(state_I1, params)

    """
    I2: NOT A1 OR NOT B1
    """
    state_I2 = A1, B1, L_A1_1, L_B1_1, N_A1_2, N_B1_2, I2
    dL_A1_1, dL_B1_1, dN_A1_2, dN_B1_2, dI2 = not_not_or2(state_I2, params)

    """
    I3: A1 OR B1
    """
    state_I3 = A1, B1, N_A1_3, N_B1_3, I3
    dN_A1_3, dN_B1_3, dI3 = yes_yes_or2(state_I3, params)

    """
    I4: NOT A0 OR NOT B0
    """
    state_I4 = A0, B0, L_A0_0, L_B0_0, N_A0_0, N_B0_0, I4
    dL_A0_0, dL_B0_0, dN_A0_0, dN_B0_0, dI4 = not_not_or2(state_I4, params)

    """
    I5: NOT I0 OR NOT I1 OR NOT A0 OR NOT B0
    """
    state_I5 = I0, I1, A0, B0, L_I0, L_I1, L_A0_1, L_B0_1, N_I0, N_I1, N_A0_1, N_B0_1, I5
    dL_I0, dL_I1, dL_A0_1, dL_B0_1, dN_I0, dN_I1, dN_A0_1, dN_B0_1, dI5 = not_not_not_not_or4(state_I5, params)

    """
    I6: NOT I2 OR NOT I3
    """
    state_I6 = I2, I3, L_I2, L_I3, N_I2, N_I3, I6
    dL_I2, dL_I3, dN_I2, dN_I3, dI6 = not_not_or2(state_I6, params)

    """
    I7: NOT I6 OR NOT I4
    """
    state_I7 = I6, I4, L_I6, L_I4, N_I6, N_I4, I7
    dL_I6, dL_I4, dN_I6, dN_I4, dI7 = not_not_or2(state_I7, params)

    """
    I8: NOT I5 OR NOT I7
    """
    state_I8 = I5, I7, L_I5, L_I7, N_I5, N_I7, I8
    dL_I5, dL_I7, dN_I5, dN_I7, dI8 = not_not_or2(state_I8, params)

    """
    S1: NOT I8
    """
    state_S1 = L_I8, S1, I8, N_I8
    dL_I8, dS1 = not_cell_wrapper(state_S1, params_not)
    dN_I8 = 0

    dA0, dA1, dB0, dB1 = 0, 0, 0, 0

    return dA0, dA1, dB0, dB1, \
            dI0, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, \
            dL_A0_0, dL_A0_1, dL_A1_0, dL_A1_1, dL_B0_0, dL_B0_1, dL_B1_0, dL_B1_1, dL_I0, dL_I1, dL_I2, dL_I3, dL_I4, dL_I5, dL_I6, dL_I7, dL_I8, \
            dN_A0_0, dN_A0_1, dN_A1_0, dN_A1_1, dN_A1_2, dN_A1_3, dN_B0_0, dN_B0_1, dN_B1_0, dN_B1_1, dN_B1_2, dN_B1_3, dN_I0, dN_I1, dN_I2, dN_I3, dN_I4, dN_I5, dN_I6, dN_I7, dN_I8,\
            dS1



"""
ALU components cells
"""

def adder_cell(state, params):
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X = params
    params_not = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x

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
    S0, S1, C = state
    

    # S0
    state_S0 = A0, B0, \
                S0_I0, S0_I1, S0_J0, \
                S0_L_A, S0_L_B, S0_L_I0, S0_L_I1, S0_L_J0, \
                S0_N_A, S0_N_B, S0_N_I0, S0_N_I1, S0_N_J0, \
                S0

    dA_S0, dB_S0, \
    dI0_S0, dI1_S0, dJ0_S0, \
    dL_A_S0, dL_B_S0, dL_I0_S0, dL_I1_S0, dL_J0_S0, \
    dN_A_S0, dN_B_S0, dN_I0_S0, dN_I1_S0, dN_J0_S0, \
    dS0 = sum0(state_S0, params)


    # S1
    state_S1 = A0, A1, B0, B1, \
                S1_I0, S1_I1, S1_I2, S1_I3, S1_I4, S1_I5, S1_I6, S1_I7, S1_I8, \
                S1_L_A0_0, S1_L_A0_1, S1_L_A1_0, S1_L_A1_1, S1_L_B0_0, S1_L_B0_1, S1_L_B1_0, S1_L_B1_1, S1_L_I0, S1_L_I1, S1_L_I2, S1_L_I3, S1_L_I4, S1_L_I5, S1_L_I6, S1_L_I7, S1_L_I8, \
                S1_N_A0_0, S1_N_A0_1, S1_N_A1_0, S1_N_A1_1, S1_N_A1_2, S1_N_A1_3, S1_N_B0_0, S1_N_B0_1, S1_N_B1_0, S1_N_B1_1, S1_N_B1_2, S1_N_B1_3, S1_N_I0, S1_N_I1, S1_N_I2, S1_N_I3, S1_N_I4, S1_N_I5, S1_N_I6, S1_N_I7, S1_N_I8,\
                S1

    dA0_S1, dA1_S1, dB0_S1, dB1_S1, \
    dI0_S1, dI1_S1, dI2_S1, dI3_S1, dI4_S1, dI5_S1, dI6_S1, dI7_S1, dI8_S1, \
    dL_A0_0_S1, dL_A0_1_S1, dL_A1_0_S1, dL_A1_1_S1, dL_B0_0_S1, dL_B0_1_S1, dL_B1_0_S1, dL_B1_1_S1, dL_I0_S1, dL_I1_S1, dL_I2_S1, dL_I3_S1, dL_I4_S1, dL_I5_S1, dL_I6_S1, dL_I7_S1, dL_I8_S1, \
    dN_A0_0_S1, dN_A0_1_S1, dN_A1_0_S1, dN_A1_1_S1, dN_A1_2_S1, dN_A1_3_S1, dN_B0_0_S1, dN_B0_1_S1, dN_B1_0_S1, dN_B1_1_S1, dN_B1_2_S1, dN_B1_3_S1, dN_I0_S1, dN_I1_S1, dN_I2_S1, dN_I3_S1, dN_I4_S1, dN_I5_S1, dN_I6_S1, dN_I7_S1, dN_I8_S1,\
    dS1 = sum1(state_S1, params)


    # C
    state_C = A0, A1, B0, B1, \
                C_I0, C_I1, C_I2, C_I3, C_I4, C_I5, \
                C_L_A0_0, C_L_A0_1, C_L_B0_0, C_L_B0_1, C_L_I0, C_L_I1, C_L_I2, C_L_I3, C_L_I4, C_L_I5, \
                C_N_A0_0, C_N_A0_1, C_N_A1_0, C_N_A1_1, C_N_B0_0, C_N_B0_1, C_N_B1_0, C_N_B1_1, C_N_I0, C_N_I1, C_N_I2, C_N_I3, C_N_I4, C_N_I5, \
                C
    
    dA0_C, dA1_C, dB0_C, dB1_C, \
    dI0_C, dI1_C, dI2_C, dI3_C, dI4_C, dI5_C, \
    dL_A0_0_C, dL_A0_1_C, dL_B0_0_C, dL_B0_1_C, dL_I0_C, dL_I1_C, dL_I2_C, dL_I3_C, dL_I4_C, dL_I5_C, \
    dN_A0_0_C, dN_A0_1_C, dN_A1_0_C, dN_A1_1_C, dN_B0_0_C, dN_B0_1_C, dN_B1_0_C, dN_B1_1_C, dN_I0_C, dN_I1_C, dN_I2_C, dN_I3_C, dN_I4_C, dN_I5_C, \
    dC = carry_out(state_C, params)



    dA0, dA1, dB0, dB1 = 0, 0, 0, 0

    return dA0, dA1, dB0, dB1, \
        dI0_S0, dI1_S0, dJ0_S0, \
        dI0_S1, dI1_S1, dI2_S1, dI3_S1, dI4_S1, dI5_S1, dI6_S1, dI7_S1, dI8_S1, \
        dI0_C, dI1_C, dI2_C, dI3_C, dI4_C, dI5_C, \
        dL_A_S0, dL_B_S0, dL_I0_S0, dL_I1_S0, dL_J0_S0, \
        dL_A0_0_S1, dL_A0_1_S1, dL_A1_0_S1, dL_A1_1_S1, dL_B0_0_S1, dL_B0_1_S1, dL_B1_0_S1, dL_B1_1_S1, dL_I0_S1, dL_I1_S1, dL_I2_S1, dL_I3_S1, dL_I4_S1, dL_I5_S1, dL_I6_S1, dL_I7_S1, dL_I8_S1, \
        dL_A0_0_C, dL_A0_1_C, dL_B0_0_C, dL_B0_1_C, dL_I0_C, dL_I1_C, dL_I2_C, dL_I3_C, dL_I4_C, dL_I5_C, \
        dN_A_S0, dN_B_S0, dN_I0_S0, dN_I1_S0, dN_J0_S0, \
        dN_A0_0_S1, dN_A0_1_S1, dN_A1_0_S1, dN_A1_1_S1, dN_A1_2_S1, dN_A1_3_S1, dN_B0_0_S1, dN_B0_1_S1, dN_B1_0_S1, dN_B1_1_S1, dN_B1_2_S1, dN_B1_3_S1, dN_I0_S1, dN_I1_S1, dN_I2_S1, dN_I3_S1, dN_I4_S1, dN_I5_S1, dN_I6_S1, dN_I7_S1, dN_I8_S1,\
        dN_A0_0_C, dN_A0_1_C, dN_A1_0_C, dN_A1_1_C, dN_B0_0_C, dN_B0_1_C, dN_B1_0_C, dN_B1_1_C, dN_I0_C, dN_I1_C, dN_I2_C, dN_I3_C, dN_I4_C, dN_I5_C, \
        dS0, dS1, dC


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
    dL_A1, dL_M0, dN_A1, dN_M0, dI0 = not_not_or2(state_I0, params)

    """
    I1: NOT M1 OR NOT A0
    """
    state_I1 = A0, M1, L_A0, L_M1, N_A0, N_M1, I1
    dL_A0, dL_M1, dN_A0, dN_M1, dI1 = not_not_or2(state_I1, params)
    

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

    # Vhodi
    A0, A1, B0, B1, \
    I0, I1, I0_0, I1_0,\
    L_A0, L_A1, L_B0, L_B1, L_I0, L_I1, L_I0_0, L_I1_0,\
    N_A0, N_A1, N_B0, N_B1, N_I0, N_I1, N_I0_0, N_I1_0,\
    S0, S1 = state

    """
    I0: A1 OR NOT B1
    """
    state_I0 = A1, B1, L_B1, N_A1, N_B1, I0
    dL_B1, dN_A1, dN_B1, dI0 = yes_not_or2(state_I0, params)


    """
    I1: NOT A1 OR B1
    """
    state_I1 = B1, A1, L_A1, N_B1, N_A1, I1
    dL_A1, dN_B1, dN_A1, dI1 = yes_not_or2(state_I1, params)


    """
    I0_0: NOT I0 OR NOT I1 OR NOT B0 OR A0
    """
    state_I0_0 = I0, I1, B0, A0, \
                L_I0, L_I1, L_B0, \
                N_I0, N_I1, N_B0, N_B1, \
                I0_0
    dL_I0, dL_I1, dL_B0, dN_I0, dN_I1, dN_B0, dN_A0, dI0_0 = not_not_not_yes_or4(state_I0_0, params)

    """
    I1_0: NOT I1 OR NOT I0 OR NOT A0 OR B0
    """
    state_I1_0 = I1, I0, A0, B0, \
                L_I1, L_I0, L_A0, \
                N_I1, N_I0, N_A0, N_B0, \
                I1_0
    dL_I1, dL_I0, dL_A0, dN_I1, dN_I0, dN_A0, dN_B0, dI1_0 = not_not_not_yes_or4(state_I1_0, params)


    """
    S0: NOT I0 OR NOT I0_0
    """
    state_S0 = I0, I0_0, L_I0, L_I0_0, N_I0, N_I0_0, S0
    dL_I0, dL_I0_0, dN_I0, dN_I0_0, dS0 = not_not_or2(state_S0, params)

    """
    S1: NOT I1 OR NOT I1_0
    """
    state_S1 = I1, I1_0, L_I1, L_I1_0, N_I1, N_I1_0, S1
    dL_I1, dL_I1_0, dN_I1, dN_I1_0, dS1 = not_not_or2(state_S1, params)


    dA0, dA1, dB0, dB1 = 0, 0, 0, 0

    return dA0, dA1, dB0, dB1, \
            dI0, dI1, dI0_0, dI1_0,\
            dL_A0, dL_A1, dL_B0, dL_B1, dL_I0, dL_I1, dL_I0_0, dL_I1_0,\
            dN_A0, dN_A1, dN_B0, dN_B1, dN_I0, dN_I1, dN_I0_0, dN_I1_0,\
            dS0, dS1



"""
ALU component models
"""

def adder_model(state, T, params):
    return np.array(adder_cell(state, params))


def shifter_model(state, T, params):
    return np.array(shifter_cell(state, params))


def comparator_model(state, T, params):
    return np.array(comparator_cell(state, params))


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

def not_not_not_yes_or4_model(state, T, params):

    dL_A, dL_B, dL_C, dN_A, dN_B, dN_C, dN_D, dOR_out = not_not_not_yes_or4(state, params)
    dA, dB, dC, dD = 0, 0, 0, 0

    return np.array([
        dA, dB, dC, dD,
        dL_A, dL_B, dL_C, 
        dN_A, dN_B, dN_C, dN_D, 
        dOR_out
    ])

def not_not_not_not_or4_model(state, T, params):
    dL_A, dL_B, dL_C, dL_D, dN_A, dN_B, dN_C, dN_D, dOR_out = not_not_not_not_or4(state, params)
    dA, dB, dC, dD = 0, 0, 0, 0

    return np.array([
        dA, dB, dC, dD,
        dL_A, dL_B, dL_C, dL_D,
        dN_A, dN_B, dN_C, dN_D, 
        dOR_out
    ])

def not_not_not_or3_model(state, T, params):
    dL_A, dL_B, dL_C, dN_A, dN_B, dN_C, dOR_out = not_not_not_or3(state, params)
    dA, dB, dC = 0, 0, 0
    return np.array([
        dA, dB, dC,
        dL_A, dL_B, dL_C, 
        dN_A, dN_B, dN_C, 
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
    return not_not_not_not_or4_model(state, T, params)