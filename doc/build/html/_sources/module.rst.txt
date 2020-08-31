.. _module:

Analyzing Data Using Modules
===========================

The analysis script package is composed of a ``backend`` folder with many individual functions within modules and a few prepared functions outside of the ``backend`` folder that puts together the individual functions in order to readily analyze the exported data (just input the exported data file and a few parameters about the experiment). 

The structure of the package is as follows:

| DLSuR
|   └── dlsmicro
|       ├── analyze_conditions.py
|       ├── analyze_replicates.py
|       ├── analyze_time_points.py
|       ├── plot_conditions.py
|       ├── plot_replicates.py
|       ├── plot_time_points.py
|       └── backend
|           ├── analysis_tools.py
|           ├── fit_funcs.py
|           ├── io.py
|           ├── plot_tools.py
|           └── utils.py
|

Prepared Functions
------------------

For full details on the functions within each module, see the section on :ref:`Analyzing Data Using Prepared Functions <analysis>`.

Modules
-------

The modules each have multiple functions within them. These functions are grouped such that each module's functions are used towards similar goals. The modules are broken into:

    1. :ref:`Data Analysis Tools <rst-analysistool>`
    2. :ref:`Fit Functions <rst-fitfunc>`
    3. :ref:`Reading Exported Data <rst-io>`
    4. :ref:`Plotting Tools <rst-plottools>`
    5. :ref:`Utilities <rst-util>`

In the following sections, the module and its role in the process of analyzing the exported data from the Malvern Zetasizer will be explained.

.. _rst-analysistool:

Data Analysis Tools
-------------------

The functions in this module are used for the analysis of the scattering function outputted by the DLS instrument. This includes all-encompassing function `full_dlsur_analysis` as well as all of the smaller functions called on by this all-encompassing function.

When using the functions in this module, one can call on `full_dlsur_analysis` to do the entire analysis. Or, one can switch up certain analysis functions. For example, the `full_dlsur_analysis` function uses the power-law analysis method for evaluating the mean-squared displacement using `msd_local_pwr_law` and `shear_modulus`. So, if one wants to see the Laplace transform method, one can write a custom script that ties together the individual functions, using `shear_modulus_laplace_transform`. 

.. _rst-fitfunc:

Fit Functions
-------------

The functions in this module are used for fitting the scattering autocorrelation function. It consists of two different types of fit: stretched exponential (`stretched_exp`) and double exponential (`expexp`). Depending on how you want to fit the scattering autocorrelation function, you can choose which fit to use.

.. _rst-io:

Reading Exported Data
---------------------

This module contains one function, which reads the data in the exported file and outputs the data in a dictionary format for easy access by future analysis functions.

The export template specified in :ref:`Exporting Data from Malvern Zetasizer Instrument <export>` is only required if one uses the function in this module.

.. _rst-plottools:

Plotting Tools
--------------

This module contains all the functions that are needed for plotting the data from the outputted Dataframes. These functions include ones that put the data into the format needed for plotting (such as `df_to_matrix` for combining multiple replicates and determining the standard deviation), ones that plot the data (such as `plot_replicates_from_df`), and ones that add features to the plot (such as `add_w_scaling`).

.. _rst-util:

Utilities
---------

In this module, there are functions that do miscellaneous work supporting the functions in other modules. For example, there is a function that calculates the Laplace transform. These functions can be leveraged by the individual user to create custom analysis functions.

Examples of Using Prepared Functions
------------------------------------

This is covered in the :ref:`Tutorial Notebook Examples <Notebook>` section.