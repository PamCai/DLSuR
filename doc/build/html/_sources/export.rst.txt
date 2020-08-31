.. _export:

Exporting Data from Malvern Zetasizer Instrument
================================================

To use the analysis script package documented here, the data collected from the Malvern Zetasizer instrument must be exported in a specified format. 

Here, the format in which the data must be exported is detailed.

Setting up the export template
------------------------------

The data, which should be composed of one longer measurement that collects the time-averaged scattering and several shorter measurements that collect the position-dependent scattering, must be exported in a specified format, which can be made into a template to be saved in the ``Export Templates`` folder.

To do so, select all rows of the data, then go to ``File`` and then ``Export``. In the first tab, be sure to save the file where you want to save it and ensure that the file is a ``.csv`` file. In the second tab, select the option to export using a template and create a new one.

As previously documented (insert reference to paper), the new template should contain the following information in exactly this order:

1. Record Number
2. Sample Name
3. Measurement Position
4. Correlation Data
5. Correlation Delay Times
6. Distribution Fit Data
7. Distribution Fit Delay Times
8. Cumulants Fit
9. Cumulants Fit Delay Times
10. Derived Count Rate
11. Measured Intercept
12. Measured Size Baseline

Additionally, it is critical to select the option to keep each piece of data separated by tabs.

Give this template a name (not Microrheology because that is a default template).

Exporting the data
------------------

Once you have saved this template and selected it as the template for exporting your rows of data, click ``OK`` and you should see a new ``.csv`` file in the folder that you saved it to.

In the future, you can cut short your steps for exporting your data to:

1. Select all rows of your data
2. Go to ``File`` and select ``Export``.
3. Select the destination for your data and ensure the file is saved as a ``.csv`` file.
4. In the second tab, select the option to save the data using a template and choose the template you made above.
5. Click ``OK``. Your data is saved!
