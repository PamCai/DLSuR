""" Module for parsing data exported from Zetasizer software"""
import pandas as pd
import numpy as np

# Columns name order for the dlsmicro_export.edf template
default_column_order = ['Record', 'Sample Name', 'Measurement Position',
                        'Correlation Data', 'Correlation Delay Times',
                        'Distribution Fit Data',
                        'Distribution Fit Delay Times',
                        'Cumulants Fit Data', 'Cumulants Fit Delay Times',
                        'Derived Count Rate', 'Measured Intercept',
                        'Measured Baseline']


def read_zetasizer_csv_to_dict(file_path, row,
                               intensities_rows=None,
                               column_order=default_column_order,
                               use_zetasizer_g1=True):
    """ Read csv file exported from the Zetasizer software to a
    dictionary containing data relevant to DLS microrheology analysis

    Parameters
    ----------
    file_path : str
                Path to the .csv file to be read
    row : int
          Row number (0-indexed) for the measurement record containing the
          correlation data
    intensities_rows : list of int, `optional`
                       List of rows (0-indexed) corresponding to the
                       scattering intensity measurements for the
                       broken-ergodicity
                       correction. If ``None``, it is assumed that the
                       intensity measurements
                       begin after the row number for the correlation data and
                       end at the last
                       row of the .csv file.
    column_order : list of str, `optional`
                   Ordered list names for the columns in the .csv file
                   (depends on your export template). If you use
                   the `dlsmicro_export.edf`, this parameter is not necessary.
                   See ``dlsmicro.io.default_column_order``
    use_zetasizer_g1 : boolean, `optional`
                       If `True`, the estimated intermediate scattering
                       function `g1`, measured
                       baseline, and measured intercept exported by the
                       Zetasizer software are used to calculate the correlation
                       function. This is useful because `g1` is exported at
                       higher numerical
                       precision than the correlation function.
                       The `g1` exported from the
                       Zetasizer is `not` the true `g1` for non-ergodic
                       samples, but this
                       option correctly inverts the formula used by
                        the Zetasizer.

    Returns
    -------
    data_dict : dictionary
                Python dictionary containing the keys below
    'time_lag' : 1d-array
                 Vector of time-lags (in microseconds) at which the
                 correlation function is measured.
    'correlation' : 1d-array
                    Vector of values of the correlation coefficient at the
                    time-lags
                    ``data_dict['time_lag']``
    'point_intensity' : float
                        Scattering intensity at the measurement position where
                        ``data_dict['correlation']`` is collected
    'ensemble_intensity' : 1d-array
    'point_position' : float
                       Meausrement position in the cuvette (in mm) where the
                       data for the
                       correlation function ``data_dict['correlation']``
                       is collected
    'ensemble_positions' : 1-d array
                           Vector of measurement positions in the cuvette
                           corresponding to
                           the scattering intensities in
                           ``data_dict['ensemble_intensity']``

    Notes
    -----
    This function will correctly parse the .csv file generated using the
    dlsmicro_export.edf Zetasizer software template. If you would like to use
    this function for parsing .csv files
    from a different user-generated Zetasizer template, 
    the following parameters `must`
    be exported. In addition,
    a list of the paramters in the order in which they occur in the .csv
    columns must be supplied. Note that
    if generating your own Zetasizer template, you should not use the "include
    headers" option. This does not work
    well for exporting the correlation data.

    'Correlation Data' : Zetasizer exports ss a string of comma separated
                         values, e.g. "1.000, 0.987, 0.921, ..."
    'Correlation Data Delay Times' : Zetasizer exports as a string of comma
                                     separated values, e.g. "0.50, 1.0, 1.5, ..."
    'Derived Count Rate'
    'Measurement Position'

    If ``use_zetsizer_g1==True``, then the export template ``must`` also
    include the following:
    'Distribution Fit Data': Zetasizer exports as a string of comma separated
     values, e.g. "0.50, 1.0, 1.5, ..."
    'Distribution Fit Delay Times': Zetasizer exports as a string of comma
                                    separated values,
                                    e.g. "0.50, 1.0, 1.5, ..."
    'Measured Intercept'
    'Measured Baseline'
    """
    # Read csv to pandas dataframe
    df = pd.read_csv(file_path, header=None, names=column_order)
    # By default, assume that scattering intensity measurements for
    # broken ergodicity correction are in second row until the end of the file
    if intensities_rows is None:
        intensities_rows = range(row + 1, len(df))

    g = np.array(
        [float(i) for i in df.iloc[row]['Correlation Data'].split(',')])
    t = np.array(
        [float(i) for i in df.iloc[row]['Correlation Delay Times'].split(',')])
    # Get the g1 correlation function data from zetasizer, which is
    # more precise than g2
    tfit = np.array(
        [float(i) for i in df.iloc[row]
         ['Distribution Fit Delay Times'].split(',')])
    g1fit = np.array(
        [float(i) for i in df.iloc[row]
         ['Distribution Fit Data'].split(',')])

    B = df.iloc[row]['Measured Baseline']
    # Get scattering intensity for the row of interest and
    # the ensemble
    Ie = df.iloc[intensities_rows]['Derived Count Rate']
    Ip = df.iloc[row]['Derived Count Rate']

    point_pos = df.iloc[row]['Measurement Position']
    epos = df.iloc[intensities_rows]['Measurement Position']
    g = np.copy(g)

    # Replace g with the data obtained from the g1 correlation
    # function where the data exists
    if use_zetasizer_g1:
        gadj = B + g1fit**2.
        tinds = [np.argmin(np.abs(t-ti)) for ti in tfit]
        g[tinds] = gadj

    data_dict = {'time_lag': t, 'correlation': g, 'point_intensity': Ip,
                 'ensemble_intensities': Ie, 'point_position': point_pos,
                 'ensemble_positions': epos}
    return data_dict
