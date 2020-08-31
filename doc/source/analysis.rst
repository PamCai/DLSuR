.. _analysis:

Analyzing Data Using Prepared Functions
=======================================

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

Modules
-------

For full details on the functions within each module, see the section on :ref:`Analyzing Data Using Modules <module>`.

Prepared Functions
------------------

For the prepared functions, where are three different kinds of data that can be analyzed and plotted (one function each for analysis and plotting).

    1. :ref:`Replicates <rst-replicates>`
    2. :ref:`Conditions <rst-conditions>`
    3. :ref:`Time Points <rst-timepoints>`

Any analysis that does not follow the format of these prepared functions can be done by piecing together a custom analysis script using the individual functions within the modules.

.. _rst-replicates:

Replicates
----------

The simplest prepared function is the one to analyze replicates of one condition. For example, if you measure the rheology 3 times of the same material, then each of those are **replicates**. So, you might use the ``analyze_replicates.py`` script to analyze multiple replicates of the same condition (or even just one replicate of one condition). 

The data should be stored in the following folder structure before running the prepared function:

| root_folder
|   ├── replicate1
|   │   └── exported.csv
|   ├── replicate2
|   │   └── exported.csv
|   └── replicate3
|       └── exported.csv
|

.. note::

    The ``analyze_replicates.py`` script requires that all replicates are named the same name (i.e. all file names are ``exported.csv``). Also, the script requires that the replicate folder names follow the above naming convention (i.e. ``replicate`` + (# of replicate)). 

Lastly, the ``plot_replicates.py`` script takes in the outputted Pandas Dataframe and plots all replicates together (saving the Dataframe must be ``True`` in ``analyze_replicates.py`` in order for ``plot_replicates.py`` to work).

.. _rst-conditions:

Conditions
----------

When there are multiple types of materials that are being measured, these types are called **conditions**. The ``analyze_conditions.py`` script is for analyzing across different material conditions.

The folder structure before running ``analyze_conditions.py`` should follow:

| root_folder
|   ├── condition1
|   │   ├── replicate1
|   │   ├── replicate2
|   │   └── replicate3
|   ├── condition2
|   │   ├── replicate1
|   │   ├── replicate2
|   │   └── replicate3
|   └── condition3
|       ├── replicate1
|       ├── replicate2
|       └── replicate3
|

.. note::

    The ``analyze_conditions.py`` script requires that all replicates are named the same name (i.e. all file names are ``exported.csv``). Also, the script requires that the replicate folder names follow the above naming convention (i.e. ``replicate`` + (# of replicate)). The ``root_folder`` and ``condition#`` folders can be anything you want, since the specific folder names must be fed into the function.

Again, the ``plot_conditions.py`` script, similarly to the ``plot_replicates.py`` script takes in the outputted Pandas Dataframe and plots all replicates and conditions together and requires the Dataframe be saved using ``analyze_conditions.py``.

.. _rst-timepoints:

Time Points
-----------

The ``analyze_time_points.py`` script is intended to facilitate measurements that are done continuously. For example, take 10 long time measurements in a row (varying in time) and one set of shorter measurements (varying in position). These types of experiments are for transient processes such as polymer gelation that occurs over the course of a couple hours to multiple hours.

The file can be inputted directly into the function without any specific file structure since all of the different time points are saved in one data file.

The ``plot_time_points.py`` script again requires that a Dataframe be saved first and plots each time point data together on one plot.

Examples of Using Prepared Functions
------------------------------------

This is covered in the :ref:`Tutorial Notebook Examples <Notebook>` section.