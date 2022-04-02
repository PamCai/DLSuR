###########################################################
# Import the functions from the DLSuR package
# You can delete lines that import functions you don't need
###########################################################
from dlsmicro.analyze_conditions import analyze_conditions
from dlsmicro.plot_conditions import plot_conditions
from dlsmicro.analyze_replicates import analyze_replicates
from dlsmicro.plot_replicates import plot_replicates
from dlsmicro.analyze_time_points import analyze_time_points
from dlsmicro.plot_time_points import plot_time_points
import numpy as np

###########################################################
# Analyzing replicates of 1 sample that you measured
###########################################################
# best to name all csv files the same
csv_name = 'exported2.csv' 
# file path to replicates
root_folder = 'example_data/replicate_example' 
# numbered files for each replicate as 'replicateX' where X = [1,2,3]
replicates = [1,2,3] 
# temperature (in Kelvin)
T = 37. + 273.15
# radius of particle (in nm)
r = 500./2.
# ergodic assumption (True if ergodic)
erg = True

# function to use to analyze the data
analyze_replicates(csv_name, root_folder, replicates, 
                T, r, erg, save_as_text=True, save_as_df=True)
# all analyzed data is saved as a Dataframe if save_as_df = True
saved_df = root_folder + '/' + 'replicate_data.pkl'
# function to plot average of analyzed data across replicates
plot_replicates(saved_df, replicates, plot_ci=True, plot_scattering=True)

##################################################################
# Analyzing replicates of 2 samples (conditions) that you measured
##################################################################
# best to name all csv files the same
csv_name = 'exported2.csv'
# file path to conditions
root_folder = 'example_data/condition_example'
# dictionary listing the name of the folders associated with each
# condition where condition 1 is the folder 'cold1', etc
condition_dir = {'condition1': 'cond1', 
				 'condition2': 'cond2'}
# numbered files for each replicate as 'replicateX' where X = [1,2,3]
# in each condition
replicate_dict = {'condition1': [1,2,3], 
                 'condition2': [1]}
# dictionary specifying colors for plotting each condition
cond_color = {'condition1': 'r', 'condition2': 'b'}
# dictionary specifying temperature of each condition
T = {'condition1': 37. + 273.15, 'condition2': 25. + 273.15}
# dictionary specifying particle size of each condition
r = {'condition1': 500./2., 'condition2': 1000./2.}
# dictionary specifying ergodicity assumption of each condition
erg = {'condition1': True, 'condition2': False}

# function to use to analyze the data
analyze_conditions(csv_name, root_folder, condition_dir, 
	               replicate_dict, T, r, erg, Laplace=True,
	               save_as_text=True, save_as_df=True,
	               plot_corr=True, plot_msd=True, plot_G=True)
# all analyzed data is saved as a Dataframe if save_as_df = True
saved_df = root_folder + '/' + 'condition_data.pkl'
# function to plot average of analyzed data across replicates for 
# all conditions
plot_conditions(saved_df, condition_dir, replicate_dict, plot_ci=True,
	            cond_color=None, plot_scattering=True, 
	            add_scaling=True, scaling_frac=[3.,4.])

##################################################################
# Analyzing time-dependent data (1 replicate of 1 sample)
##################################################################
# file path to data file
file_path = 'example_data/time_example/disposable_example.csv'
# temperature (in Kelvin)
T = 37. + 273.15
# radius of particle (in nm)
r = 30./2.
# ergodicity assumption
erg = False
# number of time points
n_points = 12
# number of positions
n_positions = 21

# function to analyze data
analyze_time_points(file_path, r, T, erg, n_points, n_positions,
                        df_save_path=None, df_file_name=None,
                        save_as_txt=True, save_as_df=True, 
                        plot_corr=False, plot_msd=False, plot_G=False)
# analyzed data saved as Dataframe if save_as_df = True
saved_df = 'example_data/time_example/time_course.pkl'
# function to plot the time-dependent rheology
plot_time_points(saved_df, n_points, plot_G_replicates=True,
					plot_MSD_replicates=True, plot_scattering=True, 
					add_scaling=True, scaling_frac=[3.,4.])