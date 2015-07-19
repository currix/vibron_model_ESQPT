def plot_components_grid(eigvec_file, min_eigstate_index, state_step, num_rows, num_cols, X_vector = 0, majorX = 50, majorY = 0.02, col="b-o"):
    ###
    '''Plot in a num_rows x num_cols grid squared components of eigstates 
    from min_eigenstate_index-th eigenvector (min value = 0) adding state_step.
    '''
    #
    import numpy as np
    from matplotlib import pyplot
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    #
    from survival_probability_basis import read_eigenstates
    #
    #
    eigstates=read_eigenstates(eigvec_file)
    #
    #
    fig,axes = pyplot.subplots(num_rows,num_cols,sharex=True,sharey=True)
    ##
#    if (not isinstance(X_vector,np.ndarray)):
#        majorXLocator   = MultipleLocator(majorX)
#        majorXFormatter = FormatStrFormatter('%d')
#    else:
#        axes.tick_params(axis='both', which='major', labelsize=12)
    ##
#    majorYLocator   = MultipleLocator(majorY)
    #
    ##pyplot.xlim(-5,240)
    #
    index_plot = min_eigstate_index
    #
    for index_row in range(0,num_rows):
        for index_col in range(0,num_cols):
            ##
            components = eigstates[:,index_plot]    ## Eigenstate's components
            ##
            axes[index_row,index_col].tick_params(axis='both', which='major', labelsize=12)
            ##
            if (not isinstance(X_vector,np.ndarray)):
                axes[index_row,index_col].plot(components**2,col) 
 #               axes[index_row,index_col].xaxis.set_major_formatter(majorXFormatter)
            else:
                axes[index_row,index_col].plot(X_vector, components**2, col, markersize = 6) 
            string_label = "$k = " + str(index_plot) + "$"
            axes[index_row,index_col].text(0.6, 0.05, string_label, fontsize=12)
            ##
            ## axes[index_row,index_col].annotate(min_eigstate_index+index_row+index_col+1, xy=(11, 0.75))
            ##
  #          axes[index_row,index_col].xaxis.set_major_locator(majorXLocator)
  #          axes[index_row,index_col].tick_params(axis='both', which='major', labelsize=9)
  #          axes[index_row,index_col].annotate(eigstate_index, xy=(150, 0.2))
            ##
            index_plot = index_plot + state_step

    pyplot.subplots_adjust(wspace=0,hspace=0)
