import numpy as np
from dlsmicro.backend import plot_tools
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rc
import matplotlib as mpl

def plot_conditions(df_path, condition_dir, replicate_dict, cond_color=None, 
					cond_label=None, plot_ci=True, plot_G_replicates=True, 
					plot_alpha_replicates=True, plot_scattering=False, 
					add_scaling=False, scaling_frac=None):

	""" Plot DLS microrheology output data for
    multiple conditions after analysis and saving as dataframe.

    Parameters
    ----------
    df_path : str
    		  Path to saved dataframe to plot
    condition_dir : dictionary
    			    Dictionary of conditions and respective folders
    cond_color : dictionary, `optional`
    			 Dictionary of conditions and respective color
    			 If `None`, colors are default colors
    cond_label : dictionary, `optional`
    			 Dictionary of conditions and respective labels
    			 If `None`, labels are taken to be folder name
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
	
	conditions = list(condition_dir.keys())
	df = pd.read_pickle(df_path)

	# set colors and labels in dictionaries
	if cond_color == None:
		cmap = plt.cm.get_cmap('hot',len(conditions)*2)
		colors = [cmap(i) for i in range(len(conditions))]
		cond_color = dict(zip(conditions, colors))
	if cond_label == None:
		cond_label = condition_dir 

	if plot_G_replicates:
		# Plot style options
		rc('axes', labelsize=24.)
		rc('xtick', labelsize=18.)
		rc('ytick', labelsize=18.)
		rc('lines', markersize=10)
		rc('lines', linewidth=3.)

		fig, ax1 = plt.subplots(1, 1)
		dash_style = (6, 2)
		for condition in conditions:
			dfi = df[df['condition'] == condition]
			plot_tools.plot_replicates_from_df(dfi, 'G1', plot_ci=plot_ci, 
											   color=cond_color[condition])
			plot_tools.plot_replicates_from_df(dfi, 'G2', plot_ci=plot_ci, 
				                               color=cond_color[condition], ls='--')

		if add_scaling:
			plot_tools.add_w_scaling(dfi[dfi['replicate'] == replicate_dict[conditions[0]][0]]['omega'][0:50],
						  scaling_frac,
						  dfi[dfi['replicate'] == replicate_dict[conditions[0]][0]]['G2'][2],
						  [2./6.,3./6.])

		# Define line objects for legend so that line color is black rather than
		# inheriting the color of one of the conditions
		handles = [None] * (len(conditions)+2)
		for con in range(len(conditions)):
			handles[con+2] = plt.Line2D((0, 1), (0, 0), color=cond_color[conditions[con]], 
										label=cond_label[conditions[con]])
		handles[0] = plt.Line2D((0, 1), (0, 0), color='k', dashes='',
		            label='$\mathregular{G^{\prime}}$')
		handles[1] = plt.Line2D((0, 1), (0, 0), color='k', ls='--',
		            label='$\mathregular{G^{\prime \prime}}$')
		plt.legend(handles=handles, frameon=False, fontsize=10, handlelength=2.0)
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
		for condition in conditions:
			dfi = df[df['condition'] == condition]
			plot_tools.plot_replicates_from_df(dfi, 'alpha', plot_ci=plot_ci, 
											   color=cond_color[condition])

		y_matrix = plot_tools.df_to_matrix(df[df['condition'] == 'condition1'], 
			                               'alpha','replicate')
		replicates = df['replicate'].values
		time = df[df['replicate'] == replicates[0]]['omega'].values
		plt.plot(time,np.full(np.shape(time),2./3.),'k--')

		# Define line objects for legend so that line color is black rather than
		# inheriting the color of one of the conditions
		handles = [None] * len(conditions)
		for con in range(len(conditions)):
			handles[con] = plt.Line2D((0, 1), (0, 0), color=cond_color[conditions[con]], 
									  label=cond_label[conditions[con]])
		plt.legend(handles=handles, frameon=False, fontsize=10, handlelength=2.0)
		plt.ylabel('$\\mathregular{\\alpha}$')
		plt.xlabel('$\\mathregular{\\omega\ (s^{-1})}$')
		plt.xscale('log')
		plt.tight_layout()
		plt.show()

	if plot_scattering:
		# Plot style options
		rc('axes', labelsize=14.)
		rc('xtick', labelsize=14.)
		rc('ytick', labelsize=14.)
		rc('lines', markersize=10)
		rc('lines', linewidth=3)

		for condition in conditions:
			dfi = df[df['condition'] == condition]
			ids = set(dfi['replicate'].values)
			for idx in ids:
				epos = dfi[dfi['replicate'] == idx]['epos'].values
				scattering = dfi[dfi['replicate'] == idx]['scattering'].values
				plot_over = np.argwhere(epos)
				plt.plot(epos[plot_over],scattering[plot_over],color=cond_color[condition],
					     linestyle='None',marker='.')

		# Define line objects for legend so that line color is black rather than
		# inheriting the color of one of the conditions
		handles = [None] * len(conditions)
		for con in range(len(conditions)):
			handles[con] = plt.Line2D((0, 1), (0, 0), color=cond_color[conditions[con]],
		                       ls='None', marker='o', label=cond_label[conditions[con]])
		plt.legend(handles=handles, frameon=False, fontsize=10, handlelength=2.0)
		plt.ylabel('Scattering Intensity')
		plt.xlabel('Position')
		plt.yscale('log')
		plt.tight_layout()
		plt.show()
