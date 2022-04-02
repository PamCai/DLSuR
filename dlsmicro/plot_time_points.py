import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from dlsmicro.backend import plot_tools
import pandas as pd

def plot_time_points(df_path, n_points, plot_G_replicates=True,
					plot_MSD_replicates=True, plot_scattering=False, 
					add_scaling=False, scaling_frac=None):

	""" Plot DLS microrheology output data for time-dependent 
	measurements after analysis and saving as dataframe.

    Parameters
    ----------
    df_path : str
    		  Path to saved dataframe to plot
    n_points : int
               Corresponds to the number of measurements taken over 
               duration of experiment
    plot_G_replicates : boolean, `optional`
    					If `True`, plot shear modulus for all time points
    plot_MSD_replicates : boolean, `optional`
    					  If `True`, plot mean-squared displacement
    					  for all time points
    plot_scattering : boolean, `optional`
    				  If `True`, plot scattering for all time points
    add_scaling : boolean, `optional`
    			  If `True`, add scaling defined by `scaling_frac` to plot
    scaling_frac : list of float, `optional`
    			   List of 2 floats, where the first number is numerator 
    			   of fraction and second is denominator of fraction
    """

	# Plot style options
	rc('axes', labelsize=20.)
	rc('xtick', labelsize=16.)
	rc('ytick', labelsize=16.)
	rc('lines', markersize=10)
	rc('lines', linewidth=1)

	# Row numbers for correlation time points (time_points)
	time_points = range(1, n_points, 1)
	# if cuvette == 'disposable':
	# 	time_points = range(1, 12, 1)
	# elif cuvette == 'quartz':
	# 	time_points = range(1, 18, 1)
	# else:
	# 	raise Exception('cuvette type entered is not valid')

	cmap = plt.cm.get_cmap('viridis',len(time_points))
	colors = [cmap(i) for i in range(len(time_points))]

	df = pd.read_pickle(df_path)

	if plot_MSD_replicates:
		for i, tp in enumerate(time_points):
			dfi = df[df['time_point'] == tp]
			time = dfi['t']
			msd = dfi['msd_smooth']
			plt.plot(time*1.e6, msd, color=colors[i], ls='-')
		plt.legend(time_points, frameon=False)
		plt.xscale('log')
		plt.yscale('log')
		plt.xlabel('$\\mathregular{\\tau\ (s)}$')
		plt.ylabel('$\\mathregular{MSD\ (nm^2)}$')
		plt.tight_layout()
		plt.show()

	if plot_G_replicates:
		for i, tp in enumerate(time_points):
			dfi = df[df['time_point'] == tp]
			freq = dfi['omega']
			G1 = dfi['G1']
			G2 = dfi['G2']
			plt.plot(freq, G1, color=colors[i], ls='-')
			plt.plot(freq, G2, color=colors[i], ls='--')
		if add_scaling:
			plot_tools.add_w_scaling(df[df['time_point'] == time_points[0]]['omega'][0:50],
						  scaling_frac,
						  df[df['time_point'] == time_points[0]]['G2'][2],
						  [2./6.,3./6.])
		plt.legend(time_points, frameon=False)
		plt.xscale('log')
		plt.yscale('log')
		plt.xlabel('$\\mathregular{\\omega\ (s^{-1})}$')
		plt.ylabel('$\\mathregular{G^*\ (Pa)}$')
		plt.tight_layout()
		plt.show()

	if plot_scattering:
		dfi = df[df['time_point'] == time_points[0]]
		epos = dfi['epos'].values
		scattering = dfi['scattering'].values
		plot_over = np.argwhere(epos)
		plt.plot(epos[plot_over],scattering[plot_over],linestyle='None',marker='.')
		plt.ylabel('Scattering Intensity')
		plt.xlabel('Position')
		plt.yscale('log')
		plt.tight_layout()
		plt.show()

