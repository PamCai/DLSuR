import numpy as np
import matplotlib.pyplot as plt

def df_to_matrix(df, quantity, replicate_identifier):

	""" Construct a matrix from a Dataframe for a given vector-valued 
	quantity in which each vector replicate is labeled by a unique 
	replicate identifier.

    Parameters
    ----------
    df : DataFrame
         Dataframe containing table of results from DLS microrheology
         analysis for a single condition
    quantity : str
          	   Name of variable to plot as defined in the Dataframe
   	replicate_identifier : str
   				 		   Name of quantity to average over

   	Returns
    -------
    M : 2-d array
        Matrix where each row is a replicate of the vector quantity.
    """

    # Construct a matrix in which each row represents a frequency sweep
	ids = set(df[replicate_identifier].values)
	M_list = []
	list_len = []
	for idx in ids:
	    xi = df[df[replicate_identifier] == idx][quantity].values
	    list_len.append(len(xi))
	    M_list.append(xi)
    # Ensure all rows are the same length
	limit = np.min(np.array(list_len))
	M_list = list(m[0:limit] for m in M_list)
	M = np.vstack(M_list)
	return M

def bootstrap_matrix_byrows(M, n_bootstrap, estimator):

	""" Gets bootstrap samples of an estimator for frequency (or time)
	sweep data from a matrix containing all vectors for a given quantity
	over all replicates.

    Parameters
    ----------
    M : 2-d array
        Matrix of all values for a given quantity over all replicates
   	n_bootstrap : int
   				  Number of points for bootstrap
    estimator : callable function
    			Function for evaluating center of distribution

   	Returns
    -------
    M_bootstrap : 2-d array
        	 	  Matrix in which each row represents a frequency sweep
    """
	n_rep = M.shape[0]
	M_bootstrap = np.zeros((n_bootstrap, M.shape[1]))
	M_sample = np.zeros(M.shape)
	for i in range(n_bootstrap):
	    inds = np.random.randint(0, n_rep, n_rep)
	    M_sample = M[inds, :]
	    M_bootstrap[i, :] = estimator(M_sample, axis=0)
	return M_bootstrap

def bootstrap_freq_sweep(df, quantity, replicate_identifier,
                         n_bootstrap, estimator=np.mean):

	""" Gets bootstrap samples of an estimator for frequency (or time)
	sweep data from a Dataframe. The Dataframe is assumed to contain a 
	vector for a given quantity in which each vector replicate is 
	labeled by a unique replicate identifier.

    Parameters
    ----------
    df : DataFrame
         Dataframe containing table of results from DLS microrheology
         analysis for a single condition
    quantity : str
          	   Name of variable to plot as defined in the Dataframe
   	replicate_identifier : str
   				 		   Name of quantity to average over
   	n_bootstrap : int
   				  Number of points for bootstrap
    estimator : callable function
    			Function for evaluating center of distribution

   	Returns
    -------
    M_bootstrap : 2-d array
        	 	  Matrix in which each row represents a frequency sweep
    """

	ids = set(df[replicate_identifier].values)
	M_list = []
	list_len = []
	for idx in ids:
	    xi = df[df[replicate_identifier] == idx][quantity].values
	    list_len.append(len(xi))
	    M_list.append(xi) 
	limit = np.min(np.array(list_len))
	M_list = list(m[0:limit] for m in M_list)
	M = np.vstack(M_list)

    # Get a matrix of bootstrapped row-wise averages given by the estimator
	M_bootstrap = bootstrap_matrix_byrows(M, n_bootstrap, estimator)

	return M_bootstrap

def bootstrap_freq_sweep_ci(df, quantity, replicate_identifier,
                            n_bootstrap, ci, estimator=np.mean):

	""" Gets bootstrap confidence interval for an estimator of
	frequency sweep (or time sweep) data. Boot strap can be either a 
	percentile bootstrap of an estimator of a studentized bootstrap
	of the mean

    Parameters
    ----------
    df : DataFrame
         Dataframe containing table of results from DLS microrheology
         analysis for a single condition
    quantity : str
          	   Name of variable to plot as defined in the Dataframe
   	replicate_identifier : str
   				 		   Name of quantity to average over
   	n_bootstrap : int
   				  Number of points for bootstrap
   	ci : float
   		 Percent of distribution included in error bars
    estimator : callable function, `optional`
    			Function for evaluating center of distribution
    
   	Returns
    -------
    ci_low : 1-d array
        	 Vector of the lower bound of the confidence interval over 
        	 entire frequency range.
    ci_high : 1-d array
			  Vector of the upper bound of the confidence interval
			  over entire frequency range.
    """

	M_bootstrap = bootstrap_freq_sweep(df, quantity, replicate_identifier,
                                       n_bootstrap, estimator=estimator)
	ci_low = np.percentile(M_bootstrap, 50.-ci/2., axis=0)
	ci_high = np.percentile(M_bootstrap, 50.+ci/2., axis=0)
	return [ci_low, ci_high]

def plot_replicates_from_df(df, my_quantity, plot_ci=True, myci=68., 
							estimator=np.mean, color='m', ls='-', 
							err_alpha=0.25, err_lw=2.5, identifier='replicate'):

    """ Plot a given quantity from the Dataframe, averaging across
    all replicates in that Dataframe.

    Parameters
    ----------
    df : DataFrame
         Dataframe containing table of results from DLS microrheology
         analysis for a single condition
    my_quantity : str
          		  Name of variable to plot as defined in the Dataframe
    plot_ci : boolean, `optional`
        	  If `True`, plot error bars
    myci : float, `optional`
    	   Percent of distribution plotted in error bars 
   	estimator : callable function, `optional`
   				Function for evaluating main plotted value
   	color : str, `optional`
   			Color of plotted line
   	ls : str, `optional`
   		 Linestyle of plotted line
   	err_alpha : float, `optional`
   				Transparency level of error bars
   	err_lw : float, `optional`
   			 Linewidth of error bar outlines
   	identifier : str, `optional`
   				 Name of quantity to average over
    """

    y_matrix = df_to_matrix(df, my_quantity, identifier)
    replicates = df['replicate'].values
    time = df[df['replicate'] == replicates[0]]['omega'].values[0:np.shape(y_matrix)[1]]
    ci = bootstrap_freq_sweep_ci(df, my_quantity, identifier, 10000,
    	                         myci, estimator=estimator)
    ci_low = ci[0]
    ci_high = ci[1]
    y_mu = estimator(y_matrix, axis=0)
    plt.plot(time, y_mu, color=color, ls=ls)
    time = np.array(time, dtype=float)
    ci_low = np.array(ci_low, dtype=float)
    ci_high = np.array(ci_high, dtype=float)
    if plot_ci:
        plt.fill_between(time, ci_low, ci_high, color=color,
                         alpha=err_alpha, linewidth=err_lw)

def add_w_scaling(omega, scaling, w_b, placement):
    
  """ Plot a given scaling on complex modulus plot to compared against the complex modulus of a sample.

    Parameters
    ----------
    omega : 1-d array
    		Vector of frequency range covered by the complex modulus
    		plotted in the plot
    scaling : list of float
          	  List of 2 floats, where the first number is numerator 
          	  of fraction and second is denominator of fraction
    w_b : float
          Value of complex modulus plotted where scaling should appear
          on the plot
    placement : list of float
    	   		First element in list is lower bound of scaling line,
    	   		second element in list is upper bound of scaling line,
    	   		where both elements are values between 0 and 1. The 
    	   		value of the first element should be less than the value
    	   		of the second element
    """

  lolim = np.int(len(omega)*placement[0])
  hilim = np.int(len(omega)*placement[1])
  omega = np.array(omega,dtype=float)
  g_scale = np.float_power(omega, scaling[0]/scaling[1])*w_b/np.float_power(omega[lolim],scaling[0]/scaling[1])
  plt.plot(omega[lolim:hilim], g_scale[lolim:hilim], ls='--', color='k',linewidth=2)
  model = np.int(0.6*(lolim+hilim))
  plt.text(omega[model],g_scale[model]*1.6,
  	     '$\omega^{%(top)s/%(bot)s}$'%{'top':np.int(scaling[0]),
  	     'bot':np.int(scaling[1])},fontsize=12)

