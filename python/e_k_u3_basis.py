def build_X_vector_U3(N_value, L_value, xi_value):
    import numpy as np
    from energy_basis_elements_I import energy_basis_I
    #
    n_vec = np.arange(L_value, N_value + 1, 2)
    #
    x_vec = []
    #
    for n in n_vec: 
        x_vec.append(energy_basis_I(N_value, n, L_value, xi_value))
    #
    return np.array(x_vec)/N_value
########
