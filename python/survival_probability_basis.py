#########################################################################
#########################################################################
def read_energy_data(filename):
    #
    import numpy as np
    #
    # Read fortran output
    fortran_output = np.loadtxt(filename)
    #
    # Select last column, eigenvalues
    eigenvalues = fortran_output[:,-1]
    #
    return eigenvalues
#########################################################################
#########################################################################
def read_eigenstates(filename):
    #
    import numpy as np
    #
    # Read fortran output
    fortran_output = np.loadtxt(filename, skiprows=1)
    #
    #
    return fortran_output
#########################################################################
#########################################################################
def survival_probability_basis_states(Eigenvalues, Eigenstates, min_index_basis, max_index_basis, t_max = 3, dim_t = 100):
    #
    '''Compute the survival probability for a given basis. See notes... (TODO complete docstring)'''
    #
    import numpy as np
    #
    time_grid = np.linspace( 0.0, t_max, dim_t)
    #
    num_curves = max_index_basis - min_index_basis + 1
    #
    result = np.zeros([num_curves,dim_t])
    #
    for curve in range(1,num_curves+1):
        #
        index_basis = curve - 1 + min_index_basis
        result[curve-1,:] = compute_surv_prob(Eigenvalues, Eigenstates, index_basis, time_grid) 
        #
        #
    return result
#########################################################################
#########################################################################
def compute_surv_prob(Energies, Eigenstates, index_basis, time_grid):
    #
    import numpy as np
    #
    surv_prob = np.zeros(len(time_grid))
    #
    # component corresponding to the index_basis state
    alpha_values = Eigenstates[index_basis-1,:] # Note that eigstates are organized in COLUMNS (file u3_2dvm_mod.f90) and that the -1 :: index_basis = 1 ==> First state
    #
    for index in range(alpha_values.shape[0]):
        surv_prob = surv_prob + alpha_values[index]**2*np.exp(Energies[index]*1j*time_grid)
    #
    return abs(surv_prob)**2 
#########################################################################
#########################################################################
def plot_surv_prob(surv_prob_output, t_max = 3, dim_t = 100, nxplot = 1, nyplot = 1):
    #
    import numpy as np
    #
    from matplotlib import pyplot
    #
    time_grid = np.linspace( 0.0, t_max, dim_t)
    #
    curve_number = surv_prob_output.shape[0]
    #
    
    pyplot.ylabel('Survival Probability')
    pyplot.xlabel('time')
    #
    for index in range(curve_number):
        pyplot.subplot(nxplot, nyplot, index)
        pyplot.plot(time_grid, surv_prob_output[index,:])
#########################################################################
#########################################################################
def energy_basis_I(N_value, n_value, l_value, xi_value):
    #
    '''Computes the diagonal matrix element <[N]n^l|H(xi)|[N]n^l> '''
    #
    return  (1-xi_value)*n_value + xi_value*(N_value*(N_value + 1) - (N_value - n_value)*(n_value + 2.0) - (N_value - n_value + 1.0)*n_value - l_value*l_value)/(N_value - 1.0)
#########################################################################
#########################################################################
def energy_basis_II(N_value, omega_value, l_value, xi_value):
    #
    '''Computes the diagonal matrix element <[N]n^l|H(xi)|[N]n^l> '''
    #
    ndiagmatel = (N_value - omega_value) * ((omega_value-l_value+2.0)*(omega_value-l_value+1.0) + (omega_value+l_value+2.0)*(omega_value+l_value+1.0)) / (2.0*(2.0*omega_value+1.0)*(2.0*omega_value+3.0)) + (N_value + omega_value + 1.0) * ((omega_value+l_value)*(omega_value+l_value-1.0) + (omega_value-l_value)*(omega_value-l_value-1.0)) / (2.0*(2.0*omega_value-1.0)*(2.0*omega_value+1.0))
    pairdiagmatel = (N_value + 1.0)*N_value - (omega_value + 1.0)*omega_value
    #
    return (1-xi_value)*ndiagmatel + xi_value*pairdiagmatel/(N_value - 1.0)
#########################################################################
#########################################################################
