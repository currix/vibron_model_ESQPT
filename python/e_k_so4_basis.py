def build_X_vector_SO4(N_value, L_value, xi_value):
    import numpy as np
    from energy_basis_elements_II import energy_basis_II
    #
    w_vec = np.arange(L_value, N_value + 1, 2)
    #
    x_vec = []
    #
    for w in w_vec: 
        x_vec.append(energy_basis_II(N_value, w, L_value, xi_value))
    #
    return np.array(x_vec)/N_value
########
