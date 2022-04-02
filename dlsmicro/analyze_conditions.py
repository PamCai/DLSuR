import numpy as np
from dlsmicro.backend import analysis_tools
from dlsmicro.backend import io
from dlsmicro.backend import utils
import pandas as pd
import matplotlib.pyplot as plt

def analyze_conditions(csv_name, root_folder, condition_dir, 
                       replicate_dict, T, r, erg, Laplace=False, 
                       df_save_path=None, df_file_name=None, 
                       save_as_text=True, save_as_df=True,
                       plot_corr=False, plot_msd=False, plot_G=False,
                       save_plots=False):

    """ Analyze files exported from Zetasizer software for multiple 
    conditions and plot data per replicate of a condition.

    Parameters
    ----------
    csv_name : str
               Name of every csv file (all the same name)
    root_folder : str
                  Name of folder in which all files to be probed are
    condition_dir : dictionary
                   Dictionary of conditions and respective folders
    replicate_dict : dictionary
                     Dictionary of replicates for each condition
    T : float or dictionary of float
        Temperature of the experiment/condition in Kelvin
    r : float or dictionary of float 
        Radius of particle in experiment/condition in nanometers
    erg : boolean or dictionary of boolean
          Ergodicity in experiment/condition
    Laplace : boolean, `optional`
              If `True`, use direct Laplace transform to find
              shear modulus. This is useful because it can
              smooth the data to noise.
    df_save_path : str, `optional` 
                   Path to Dataframe to be saved containing results
                   from DLS microrheology analysis.
                   If `None`, path will be set to root_folder.
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
    save_plots : boolean, `optional`
             If `True`, saves plots of correlation function, MSD, G
    """

    conditions = list(condition_dir.keys())
    if type(T) is dict:
        if T.keys() != condition_dir.keys():
            raise Exception('Keys for temp dictionary do not match conditions')
        T_dict = T
    else:
        T_dict = dict((k, T) for k, v in condition_dir.items())
    if type(r) is dict:
        if r.keys() != condition_dir.keys():
            raise Exception('Keys for radius dictionary do not match conditions')
        r_dict = r
    else:
        r_dict = dict((k, r) for k, v in condition_dir.items())
    if type(erg) is dict:
        if erg.keys() != condition_dir.keys():
            raise Exception('Keys for ergodicity dictionary do not match conditions')
        erg_dict = erg
    else:
        erg_dict = dict((k, erg) for k, v in condition_dir.items())
    if df_save_path == None:
        df_save_path = root_folder
    if df_file_name == None:
        df_file_name = 'condition_data.pkl'

    # Define scattering vector parameters about solvent
    n = 1.333 # index of refraction of water (default)
    theta = 173.*np.pi/180. 
    lam = 633.

    ####################################################
    # Analyze data
    # Don't edit this unless you know what you're doing
    ####################################################

    # Create pandas dataframe to organize replicate data
    df = pd.DataFrame(columns=['t', 'msd_smooth', 'alpha', 'omega', 'G1', 'G2',
                               'replicate', 'condition', 'id', 'scattering', 'epos'])
    idx = 0
    for condition in conditions:
        for replicate in replicate_dict[condition]:
            file_path = '%s/%s/replicate%s/%s' % (root_folder,
                                                  condition_dir[condition],
                                                  replicate, csv_name)

            # Read the data and truncate over the desired time windows
            data_dict = io.read_zetasizer_csv_to_dict(file_path, 0)
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
            T = T_dict[condition]
            r = r_dict[condition]
            ergodic = erg_dict[condition]
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
            dlsmicro_df['replicate'] = [replicate]*len(dlsmicro_df['t'])
            dlsmicro_df['scattering'] = scattering
            dlsmicro_df['epos'] = positions
            dlsmicro_df['condition'] = condition
            dlsmicro_df['id'] = idx

            idx += 1
            df = pd.concat((df,dlsmicro_df), axis=0, sort=True)

            ###############################################
            # Plot the analyzed data for this replicate
            ###############################################

            if plot_corr:
                plt.plot(t, g, '-r')
                plt.xscale('log')
                plt.xlabel('$\\mathregular{time\ (\\mu s)}$')
                plt.ylabel('g')
                plt.tight_layout()
                if save_plots:
                    save_path = '%s/%s/replicate%s/corr' % (root_folder,
                                                            condition_dir[condition],
                                                            replicate)
                    plt.savefig(save_path)
                plt.show()
            if plot_msd:
                plt.plot(dlsmicro_df['t'], dlsmicro_df['msd_smooth'], '-r')
                plt.xscale('log')
                plt.yscale('log')
                plt.xlabel('$\\mathregular{time\ (\\mu s)}$')
                plt.ylabel('MSD')
                plt.tight_layout()
                if save_plots:
                    save_path = '%s/%s/replicate%s/msd' % (root_folder,
                                                            condition_dir[condition],
                                                            replicate)
                    plt.savefig(save_path)
                plt.show()
            if plot_G:
                plt.plot(dlsmicro_df['omega'], dlsmicro_df['G1'], '-r')
                plt.plot(dlsmicro_df['omega'], dlsmicro_df['G2'], '--r')
                plt.ylabel('$\\mathregular{G^*\ (Pa)}$')
                plt.xlabel('$\\mathregular{\\omega\ (s^{-1})}$')
                plt.legend(['$\\mathregular{G^{\\prime}}$',
                            '$\\mathregular{G^{\\prime \\prime}}$'], frameon=False,)
                plt.xscale('log')
                plt.yscale('log')
                plt.tight_layout()
                if save_plots:
                    save_path = '%s/%s/replicate%s/G' % (root_folder,
                                                            condition_dir[condition],
                                                            replicate)
                    plt.savefig(save_path)
                plt.show()

            ##############################################
            # Save analysis results
            ##############################################
            if save_as_text:
                save_path = '%s/%s/replicate%s/' % (root_folder,
                                                    condition_dir[condition],
                                                    replicate)
                for i in dlsmicro_df.columns[:-2]:
                    np.savetxt('%s/%s' % (save_path,i), dlsmicro_df[i].values)

    #################################################
    # Save the pandas dataframe
    #################################################
    if save_as_df:
        save_path = df_save_path + '/' + df_file_name
        df.to_pickle(save_path)
