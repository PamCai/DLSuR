import numpy as np
from dlsmicro.backend import analysis_tools
from dlsmicro.backend import io
from dlsmicro.backend import utils
import matplotlib.pyplot as plt
import pandas as pd

def analyze_time_points(file_path, T, r, ergodic, n_points, n_positions,
                        Laplace=False, df_save_path=None, df_file_name=None,
                        save_as_txt=True, save_as_df=True, 
                        plot_corr=False, plot_msd=False, plot_G=False):

    """ Analyze files exported from Zetasizer software for time-
    dependent measurements and plot data per time point.

    Parameters
    ----------
    file_path : str
               Path to exported csv file
    T : float
        Temperature of the experiment in Kelvin
    r : float  
        Radius of particle in experiment in nanometers
    ergodic : boolean
              Ergodicity in experiment
    n_points : int
               Corresponds to the number of measurements taken over 
               duration of experiment
    n_positions : int
                  Corresponds to the number of measurements taken
                  over different positions at end of experiment
    Laplace : boolean, `optional`
              If `True`, use direct Laplace transform to find
              shear modulus. This is useful because it can
              smooth the data to noise.
    df_save_path : str, `optional` 
                   Path to Dataframe to be saved containing results
                   from DLS microrheology analysis.
                   If `None`, path will be set to file_path.
    df_file_name : str, `optional`
                   Name of Dataframe to be saved containing results
                   from DLS microrheology analysis.
                   If None, `condition_data.pkl` is default name
    save_as_text : boolean, `optional`
                   If `True`, save each element of Dataframe separately
                   as a text file
    save_as_df : boolean, `optional`
                 If `True`, save the Dataframe
    plot_corr : boolean, `optional`
                If `True`, show plot of the correlation function of
                each replicate
    plot_msd : boolean, `optional`
               If `True`, show plot of the mean-squared displacement 
               of each replicate
    plot_G : boolean, `optional`
             If `True`, show plot of the shear modulus of each replicate
    """

    if df_save_path == None:
        files = file_path.split('/')
        file_name = files[-1]
        df_save_path = file_path[:-(len(file_name))]
    if df_file_name == None:
        df_file_name = 'time_course.pkl'

    # Row numbers for scattering intensity sweep (int_rcds)
    # Row numbers for correlation time points (time_points)
    int_rcds = range(n_points + 1, n_points + n_positions)
    time_points = range(1, n_points, 1)
    # if cuvette == 'disposable':
    #     int_rcds = range(13, 33)
    #     time_points = range(1, 12, 1)
    # elif cuvette == 'quartz':
    #     int_rcds = range(19, 51)
    #     time_points = range(1, 18, 1)
    # else:
    #     raise Exception('cuvette type entered is not valid')

    # Define scattering vector parameters about solvent
    n = 1.333 # index of refraction of water (default)
    theta = 173.*np.pi/180. 
    lam = 633.

    # Create pandas dataframe to organize time point data
    df = pd.DataFrame(columns=['time', 'MSD', 'alpha', 'omega', 'G1', 'G2',
                               'time_point', 'scattering','epos'])

    for tp in time_points:
        save_suffix = 'time_point_%s.txt' % tp
        ####################################################
        # Analyze data
        # Don't edit this unless you know what you're doing
        ####################################################

        # Read the data and truncate over the desired time windows
        data_dict = io.read_zetasizer_csv_to_dict(file_path, tp, intensities_rows=int_rcds)
        [t, g, I, Ie, point_pos, epos] = list(data_dict.values())

        # Figure out where data is no longer trustworthy (correlation function goes to zero)
        tinds = [3, len(g)]
        tinds[1] = np.argmax(g<0.05)
        if tinds[1] == 0:
            tinds[1] = len(g) - 1
        g = g[tinds[0]:tinds[1]]
        t = t[tinds[0]:tinds[1]]

        # Get a DLS microrheology object that contains all raw data
        # and analyzed results. 
        q = analysis_tools.calc_q(n, theta, lam)
        dlsmicro_df = analysis_tools.full_dlsur_analysis(t, g, ergodic, r, T, q,
                                                         I, Ie)

        # Store the scattering vs. position data
        scattering = np.zeros(len(dlsmicro_df['t']))
        positions = np.zeros_like(scattering)
        scattering[0:len(Ie)] = Ie
        positions[0:len(epos)] = epos

        # Laplace transformed modulus
        if Laplace:
            [omega_L, G1_L, G2_L] = analysis_tools.shear_modulus_laplace_transform(t, 
                                                   dlsmicro_df['msd_smooth'], r, T)
            dlsmicro_df['G1'], dlsmicro_df['G2'] = utils.laplace_merge(dlsmicro_df['t'], 
                                                                     dlsmicro_df['G1'], 
                                                                     dlsmicro_df['G2'],
                                                                     G1_L,G2_L)

        # Construct pandas dataframe for this replicate and append it to the
        # master dataframe
        dlsmicro_df['time_point'] = [tp]*len(dlsmicro_df['t'])
        dlsmicro_df['scattering'] = scattering
        dlsmicro_df['epos'] = positions

        df = pd.concat((df,dlsmicro_df), axis=0, sort=True)

        ###############################################
        # Plot the analyzed data
        ###############################################

        if plot_corr:
            plt.plot(t, g, '-r')
            plt.xscale('log')
            plt.xlabel('$\\mathregular{time\ (\\mu s)}$')
            plt.ylabel('g')
            plt.show()
        if plot_msd:
            plt.plot(dlsmicro_df['t'], dlsmicro_df['msd_smooth'], '-r')
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel('$\\mathregular{time\ (\\mu s)}$')
            plt.ylabel('MSD')
            plt.show()
        if plot_G:
            plt.plot(dlsmicro_df['omega'], dlsmicro_df['G1'], '-r', lw=2.0)
            plt.plot(dlsmicro_df['omega'], dlsmicro_df['G2'], '--r', lw=2.0)
            plt.ylabel('$\\mathregular{G^*\ (Pa)}$')
            plt.xlabel('$\\mathregular{\\omega\ (s^{-1})}$')
            plt.legend(['$\\mathregular{G^{\\prime}}$',
                        '$\\mathregular{G^{\\prime \\prime}}$'], frameon=False,)
            plt.xscale('log')
            plt.yscale('log')
            plt.show()

        ##############################################
        # Save analysis results
        ##############################################

        if save_as_txt:
            for i in dlsmicro_df.columns:
                np.savetxt('%s/%s' % (df_save_path,i), dlsmicro_df[i].values)

    #################################################
    # Save the pandas dataframe
    #################################################
    if save_as_df:
        save_path = df_save_path + '/' + df_file_name
        df.to_pickle(save_path)
