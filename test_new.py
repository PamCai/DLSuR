from dlsmicro.analyze_conditions import analyze_conditions
from dlsmicro.plot_conditions import plot_conditions
from dlsmicro.analyze_replicates import analyze_replicates
from dlsmicro.plot_replicates import plot_replicates
from dlsmicro.analyze_time_points import analyze_time_points
from dlsmicro.plot_time_points import plot_time_points
import numpy as np

csv_name = 'exported2.csv'
root_folder = 'condition_example'
condition_dir = {'condition1': 'condition1', 
				 'condition2': 'condition2'}
replicate_dict = {'condition1': [1,2,3], 
                 'condition2': [1]}
cond_color = {'condition1': 'r', 'condition2': 'g'}

T = 37. + 273.15
r = 500./2.
erg = True

# analyze_conditions(csv_name, root_folder, condition_dir, 
# 	               replicate_dict, T, r, erg, Laplace=True,
# 	               save_as_text=True, save_as_df=False,
# 	               plot_corr=True, plot_msd=True, plot_G=True)
# plot_conditions('condition_example/condition_data.pkl',  
# 	            condition_dir, replicate_dict, plot_ci=True,
# 	            cond_color=None, plot_scattering=True, 
# 	            add_scaling=True, scaling_frac=[3.,4.])

# root_folder = 'replicate_example'
# analyze_replicates(csv_name, root_folder, replicate_dict['condition1'], 
#                 T, r, erg, save_as_text=True, 
#                 save_as_df=True)
# plot_replicates('replicate_example/replicate_data.pkl',
# 				replicate_dict['condition1'], plot_scattering=True)

file_path = 'disposable_time_example/disposable_example.csv'
r = 30./2.
ergodic = False
cuvette = 'disposable'

analyze_time_points(file_path, r, T, ergodic, cuvette=cuvette, 
                        df_save_path=None, df_file_name=None,
                        save_as_txt=True, save_as_df=True, 
                        plot_corr=False, 
                        plot_msd=False, plot_G=False)
plot_time_points('disposable_time_example/time_course.pkl', cuvette='disposable', plot_G_replicates=True,
					plot_MSD_replicates=True, plot_scattering=True, 
					add_scaling=True, scaling_frac=[3.,4.])