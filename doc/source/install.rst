.. _install:

Installation and Implementation of DLSuR Python Package
===============================================================

Installing the ``DLSuR`` Python package is very simple. It is recommended to do this within an Anaconda environment. 

A short video tutorial of creating an environment, installing DLSuR in the environment, and implementing the package to analyze data can be found below. If you would like to read the instructions instead, keeping scrolling down.

.. raw:: html

   <iframe width="560" height="315" src="https://www.youtube.com/embed/kbmUlsj4Esk" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

Instructions for each step are linked here for easy access.

    1. :ref:`Creating an Anaconda Environment <rst-create>`
    2. :ref:`Installing DLSuR Python Package <rst-install>`
    3. :ref:`Importing and Implementing DLSuR <rst-run>`

Any analysis that does not follow the format of these prepared functions can be done by piecing together a custom analysis script using the individual functions within the modules.

.. _rst-create:

Creating an Anaconda Environment
--------------------------------

To create an Anaconda environment, simple enter in the command line the following code.

.. code-block:: bash

    $ conda create --name example

You can replace ``example`` with any name you choose to name your environment.

.. _rst-install:

Installing DLSuR Python Package
-------------------------------

Before installing ``DLSuR``, enter your environment by entering into your command line the following code.

.. code-block:: bash

    $ conda activate example

Again, ``example`` should be replaced with the name you have chosen for your created environment.

Next, we need to install a preliminary package first so that we have the Python package ``pip`` installed. For this, we will install the Python package ``scipy``.

.. code-block:: bash

    $ conda install scipy

Click through the prompts (i.e. it will ask whether you want to proceed, type ``y``). Once ``scipy`` has been installed, you can double check that ``pip`` has also been installed by typing the following into your command line:

.. code-block:: bash

    $ conda list

The package ``pip`` should show up in the list that appears.

Now, we install ``DLSuR`` by entering the following code.

.. code-block:: bash

    $ pip install DLSuR

This can take some time. Please be patient. Once it is finished, congratulations! You now have ``DLSuR`` installed!

.. _rst-run:

Importing and Implementing DLSuR
--------------------------------

To actually use the functions within the ``DLSuR`` package, you need to import them into your code file. This can be done by adding the following line at the top of your code file::

    import dlsmicro

If you want to import a specific function (i.e. ``analyze_conditions``), then you can use the following syntax::

    from dlsmicro import analyze_conditions

If the above line does now allow you to find the function, try the following code::

    from dlsmicro.analyze_conditions import analyze_conditions

To use this function, define all of the inputs required, including file names, replicates list, temperature, particle radius. See below for an example of the inputs being defined::

    csv_name = 'exported2.csv'
    root_folder = 'example_data/condition_example'
    condition_dir = {'condition1': 'cond1', 'condition2': 'cond2'}
    replicate_dict = {'condition1': [1,2,3], 'condition2': [1]}
    cond_color = {'condition1': 'r', 'condition2': 'b'}
    T = {'condition1': 37. + 273.15, 'condition2': 25. + 273.15}
    r = {'condition1': 500./2., 'condition2': 1000./2.}
    erg = {'condition1': True, 'condition2': False}


The function ``analyze_conditions`` can now be used to analyze the data defined above::

    analyze_conditions(csv_name, root_folder, condition_dir, 
	               replicate_dict, T, r, erg, Laplace=True,
	               save_as_text=True, save_as_df=True,
	               plot_corr=True, plot_msd=True, plot_G=True)

The analyzed data is saved as a Pandas Dataframe in the ``root_folder`` and titled ``condition_data.pkl`` unless you have defined a different name using the input parameter ``df_file_name``. This Dataframe is the data input for the function ``plot_conditions`` that allows you to see all conditions plotted together as the average of the replicates of each condition::

    saved_df = root_folder + '/' + 'condition_data.pkl'

    plot_conditions(saved_df, condition_dir, replicate_dict, plot_ci=True,
	            cond_color=None, plot_scattering=True, 
	            add_scaling=True, scaling_frac=[3.,4.])

Please visit `this package's Github page <https://github.com/PamCai/DLSuR>`_ to find the above code blocks in the file ``test_new.py`` and the example data referenced here to see how the file structure is set up.