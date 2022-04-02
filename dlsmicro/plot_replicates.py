import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from dlsmicro.backend import plot_tools
from matplotlib import rc
import matplotlib as mpl

def plot_replicates(df_path, replicates, replic_color=None,
					plot_ci=True, plot_G_replicates=True,
					plot_alpha_replicates=True, plot_scattering=False,
					add_scaling=False, scaling_frac=None):

	""" Plot DLS microrheology output data for
    multiple replicates after analysis and saving as dataframe.

    Parameters
    ----------
    df_path : str
    		  Path to saved dataframe to plot
    replicates : list of ints
                 List of ints corresponding to the replicate number
    replic_color : dictionary, `optional`
    			 Dictionary of replicates and respective color
    			 If `None`, colors are default colors. This 
    			 applies to the scattering plot. Otherwise, color
    			 plotted for G and alpha is first color in dictionary
    plot_ci : boolean, `optional`
    		  If `True`, plot error bars around G and alpha
    plot_G_replicates : boolean, `optional`
    					If `True`, plot G for all conditions averaged
    					over all replicates
    plot_alpha_replicates : boolean, `optional`
    						If `True`, plot alpha for all conditions
    						averaged over all replicates
    plot_scattering : boolean, `optional`
    				  If `True`, plot scattering for all conditions
    add_scaling : boolean, `optional`
    			  If `True`, add scaling defined by `scaling_frac` to plot
    scaling_frac : list of float, `optional`
    			   List of 2 floats, where the first number is numerator 
    			   of fraction and second is denominator of fraction
    """

	df = pd.read_pickle(df_path)

	# set colors and labels in dictionaries
	if replic_color == None:
		cmap = plt.cm.get_cmap('hot',len(replicates)*2)
		colors = [cmap(i) for i in range(len(replicates))]
		replic_color = dict(zip(replicates, colors))

	if plot_G_replicates:
		# Plot style options
		rc('axes', labelsize=24.)
		rc('xtick', labelsize=18.)
		rc('ytick', labelsize=18.)
		rc('lines', markersize=10)
		rc('lines', linewidth=3)

		fig, ax1 = plt.subplots(1, 1)
		plot_tools.plot_replicates_from_df(df, 'G1', plot_ci=plot_ci, 
										   color=replic_color[replicates[0]])
		plot_tools.plot_replicates_from_df(df, 'G2', plot_ci=plot_ci, 
									  	   color=replic_color[replicates[0]], ls='--')

		if add_scaling:
			plot_tools.add_w_scaling(df[df['replicate'] == replicates[0]]['omega'][0:50],
				          			 scaling_frac,
				          		     df[df['replicate'] == replicates[0]]['G2'][2],
				          			 [2./6.,3./6.])
		handles = [None] * (2)
		handles[0] = plt.Line2D((0, 1), (0, 0), color=replic_color[replicates[0]], dashes='',
								label='$\mathregular{G^{\prime}}$')
		handles[1] = plt.Line2D((0, 1), (0, 0), color=replic_color[replicates[0]], ls='--',
								label='$\mathregular{G^{\prime \prime}}$')
		plt.legend(handles=handles, frameon=False, fontsize=16)
		plt.ylabel('$\\mathregular{G^*\ (Pa)}$')
		plt.xlabel('$\\mathregular{\\omega\ (s^{-1})}$')
		plt.xscale('log')
		plt.yscale('log')

		# Force plot to show ticks at every decade and minor ticks
		minor_tick_marks = [0.1*x for x in range(1, 10)]
		locmax = mpl.ticker.LogLocator(base=10, numticks=12)
		locmin = mpl.ticker.LogLocator(base=10, numticks=12, subs=minor_tick_marks)
		ax1.xaxis.set_major_locator(locmax)
		ax1.xaxis.set_minor_locator(locmin)
		ax1.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
		plt.tight_layout()
		plt.show()

	if plot_alpha_replicates:
		# Plot style options
		rc('axes', labelsize=24.)
		rc('xtick', labelsize=18.)
		rc('ytick', labelsize=18.)
		rc('lines', markersize=10)
		rc('lines', linewidth=3)

		fig, ax1 = plt.subplots(1, 1)
		plot_tools.plot_replicates_from_df(df, 'alpha', plot_ci=plot_ci, 
			                               color=replic_color[replicates[0]])
		plt.ylabel('$\\mathregular{\\alpha}$')
		plt.xlabel('$\\mathregular{\\omega\ (s^{-1})}$')
		plt.xscale('log')
		plt.tight_layout()
		plt.show()

	if plot_scattering:
		# Plot style options
		rc('axes', labelsize=18.)
		rc('xtick', labelsize=14.)
		rc('ytick', labelsize=14.)
		rc('lines', markersize=10)
		rc('lines', linewidth=3)

		fig, ax1 = plt.subplots(1, 1)

		for replicate in replicates:
			epos = df[df['replicate'] == replicate]['epos'].values
			scattering = df[df['replicate'] == replicate]['scattering'].values
			plot_over = np.argwhere(epos)
			plt.plot(epos[plot_over],scattering[plot_over],
				     color=replic_color[replicate],linestyle='None',marker='.')

		# Define line objects for legend so that line color is black rather than
		# inheriting the color of one of the conditions
		handles = [None] * len(replicates)
		for replicate in range(len(replicates)):
			handles[replicate] = plt.Line2D((0, 1), (0, 0), color=replic_color[replicates[replicate]],
											ls='None', marker='o', label=replicates[replicate])
		plt.legend(handles=handles, frameon=False, fontsize=10, handlelength=2.0)
		plt.ylabel('Scattering Intensity')
		plt.xlabel('Position')
		plt.yscale('log')
		plt.tight_layout()
		plt.show()
