""" Module for processing DLS correlation data into rheological properties"""

import numpy as np
import dlsmicro.backend.utils as utils
from scipy import special
import dlsmicro.backend.fit_funcs as fit_funcs
import pandas as pd


def find_g0(t, corr, func=fit_funcs.stretched_exp, t0=2.0,
            tmaxs=np.arange(40., 130., 10), p0=None):
    """ Estimate the intercept of the correlation function at t = 0

    Parameters
    ----------
    t  : 1d-array
         Vector of time-lags at which the correlation coefficient is computed
    corr : 1d-array
         Correlation coefficient (which spans possible from 0.0 to 1.0).
         Note that this is equal to `g2 - 1`, where `g2` is the intensity
         autocorrelation function.
    func : callable function
           Model to fit the correlation against. Must be of the form
           `f(t, p1, p2, ..., pM)` where `t` is the independent variable and
           `p1, p2, ..., pM` is a set of M parameters to fit.
    t0: float, `optional`
        Minimum time-lag to use in the fit of the correlation coefficient
        for estimating the intercept.
    tmaxs : 1-d array or list of floats, `optional`
           Vector of possible maximum time-lags to use in the fitting interval
           of the correlation coefficient. The optimal tmax will be selected
           based on minimization of the cross-validation error of the fit.
    p0 : 1-d array, `optional`
         Initial guesses for the parameters to ``func`` for fitting


    Returns
    -------
    g0 : float
         Estimate of the intercept of ``corr`` at `t=0`
    twindow_min : List of floats
                  Optimized interval in time over which the fit to ``corr``
                  was performed to obtain g0. This is a list
                  of the form [t0, tmax].
    pmin : 1d-array
           Optimal values for the fitting parameters to ``func`` used to
           estimate ``g0``
    """

    # Construct list of fitting windows
    twindows = [[t0, tmax] for tmax in tmaxs]
    # Find twindow for minimum CV error
    [twindow_min, pmin, CV_min] = utils.minimize_cv_error(t, corr, twindows,
                                                          func, p0)
    g0 = func(0.0, *pmin)

    return [g0, twindow_min, pmin]


def calc_g1(t, corr, ergodic, g0=None, Ip=None, Ie=None, eps=None):
    """Compute the intermediate scattering function from the correlation function.

    This function transforms the correlation functio to the intermediate 
    scattering function, which is necessary for computing the particle MSD.
    The function enables corrections for broken-ergodicty and estimation
    of the zero-time intercept of the correlation function.

    Parameters
    ----------
    t  : 1d-array
         Vector of time-lags at which the correlation function is computed
    corr : 1d-array
         Correlation coefficient (which spans possible from 0.0 to 1.0).
         Note that this is equal to `g2 - 1`, where `g2` is the intensity
         autocorrelation function.
    ergodic: boolean
             If False, corrections are made for non-ergodic systems
             based on ``Ip`` and ``Ie``
    Ip : float
        Scattering intensity at the position
        where correlation function is collected.
        (Only used if ``ergodic=False``)
    Ie : 1d-array
        Vector of scattering intensities at different positions in the cuvette
        (Only used if ``ergodic=False``)
    g0 : float, `optional`
         Estimate of the intercept of ``g2 - 1`` at time 0. If not provided,
         the intercept will be estimated automatically based on a stretched
         exponential fit.

    Returns
    -------
    g1 : 1d-array
        The intermediate scattering function corresponding to
        the time-lags ``t``

    Notes
    -----
    The correlation function exported by the Zetasizer software is equal to
    ``g2 - 1``
    """

    g2 = corr + 1.

    # If no intercept is provided, estimate it using a stretched exponential
    # function with default parameters

    if g0 is None:
        # Define guesses for the stretched exponential function fitting
        a0 = 1.0e-2
        beta0 = 1.0
        p0 = [corr[1], a0, beta0]
        # Fit the correlation and get the intercept
        [g0, twindow_min, pmin] = find_g0(t, corr, func=fit_funcs.stretched_exp, t0=2.0,
                     tmaxs=np.arange(40., 130., 10.), p0=p0)
    
    # If any of the g values are greater than g0, allow
    # g2 to be replaced with the fit to avoid negative
    # MSD values
    if np.any(corr > g0):
        print('Negative values encountered in MSD...')
        print('Replacing early time data with gfit')
    tmin = np.argmin(np.abs(t-twindow_min[0]))
    tmax = np.argmin(np.abs(t-twindow_min[1]))
    gfit = fit_funcs.stretched_exp(t, *pmin)
    g2[tmin:tmax] = gfit[tmin:tmax] + 1.

    if ergodic:
        g1 = np.sqrt((g2-1.)/g0)
    else:
        Ie_avg = np.average(Ie)
        # calculate the ratio of ensemble to time averaged
        # intensities
        Y = Ie_avg/Ip
        # calculate the g1 correlation function
        if eps is None:
            g1 = (Y-1.)/Y+np.sqrt(g2-g0)/Y
        else:
            g1 = (1.-(1.-eps)/Y +
                  (1-eps)*np.sqrt(1. +
                                  (g2-g0-1.)/(1-eps)**2.)/Y)
    return g1


def calc_q(n, theta, lam):
    """ Calculate the scattering wave-vector.

    This function returns the scattering wave-vector for a DLS experiment
    based on the scattering geometry and the wave-length of the laser.

    Parameters
    ----------
    n : float
        Index of refraction of the material (1.333 for water)
    theta : float
            Scattering angle in radians
    lam : float
          wave-length of laser (units used here determine the
          units of ``q``)

    Returns
    -------
    q : float
        Scattering wave-vector in inverse units of ``1/lam``
    """
    q = 4.*np.pi*n*np.sin(theta/2.)/(lam)
    return q


def msd_local_pwr_law(t, g1, q, bw=0.1, replace_neg=True):
    """ Calculate the local power-law scaling of the MSD and the
        smoothed MSD by locally-weighted logarithmic linear regression

    Parameters
    ----------
    t : 1d-array
        Vector of time-lags
    g1 : 1d-array
         Vector containing the intermediate scattering function at time-lags
         given by ``t``
    q : float
        Scattering vector in units of 1/nm
    bw : float, `optional`
           Bandwith smoothing parameter for locally-weighted regression.
           Reasonable values are typically between 0.05 and 0.1


    Returns
    -------
    msd_smooth: 1d-array
                Vector of smoothed MSD values from the local regression
                corresponding to the time lags in ``t`` (in units of nm^2).
    alpha: 1d-array
           Vector of local power-law scaling exponents of the MSD
           corresponding to the time lags in ``t``.
    """
    msd = -6*np.log(g1)/(q**2.)

    # Remove data points with 0, negative, or infinite MSD
    if replace_neg:
        if any(msd < 0):
            print('negative msd values enountered,'
                  'replacing with interpolation')
            neg_inds = msd < 0
            t_pos = t[msd > 0]
            msd_pos = msd[msd > 0]
            t_neg = t[neg_inds]
            msd_interp = np.interp(t_neg, t_pos, msd_pos)
            msd[neg_inds] = msd_interp

    # Perform a locally-weighted logarithmic linear regression
    [Theta, log_msd_smooth] = utils.loess(np.log(t), np.log(msd),
                                          degree=1, alpha=bw)
    msd_smooth = np.exp(log_msd_smooth)
    # Define the local power law scaling exponent based
    # on the logarithmic slopes
    alpha = Theta[1, :]
    return [msd_smooth, alpha]


def calc_msd_raw(t, g1, q, replace_neg=True):
    """ Calculate the MSD from the intermediate scattering function.

    Parameters
    ----------
    t : 1d-array
        Vector of time-lags
    g1 : 1d-array
         Vector containing the intermediate scattering function at time-lags
         ``t``
    q : float
        Scattering vector in units of 1/nm
    replace_neg : boolean, `optional`
                  If ``True``, replace negative values of the MSD by linear
                  interpolation

    Returns
    -------
    msd : 1d-array
          Vector of mean-squared-displacements at time-lags ``t``
          (in units of nm^2)
    """
    msd = 6*np.log(g1)/(q**2.)

    # Remove data points with 0, negative, or infinite MSD
    if replace_neg:
        if any(msd < 0):
            print('negative msd values enountered,'
                  'replacing with interpolation')
            neg_inds = msd < 0
            t_pos = t[msd > 0]
            msd_pos = msd[msd > 0]
            t_neg = t[neg_inds]
            msd_interp = np.interp(t_neg, t_pos, msd_pos)
            msd[neg_inds] = msd_interp

    return msd


def shear_modulus(t, msd, alpha, r, T):
    """ Calculate Frequency-dependent shear modulus by power-law analysis 

    Calculate the frequency-dependent complex shear modulus using results
    from a power-law analysis of the mean-squared-displacement. 

    Parameters
    ----------
    t : 1d-array
        Time-lags in units of microseconds
    msd : 1d-array
          Mean-squared displacements at the time-lags
          ``t`` in units of nm^2
    alpha : 1d-array of local power-law scaling exponents
    r : float
        Radius of the probe particles in nanometers
    T : float
        Temperature in Kelvin

    Returns
    -------
    omega : 1d-array
            Vector of angular frequencies in units of 1/s
    G1 : 1d-array
         Storage modulus at angular frequencies ``omega`` in
         units of Pa
    G2 : 1d-array
         Loss modulus at angular frequencies ``omega`` in
         units of Pa

    Notes
    -----
    This method produces more accurate results in the upper and
    lower frequency extremes than a direct numerical Laplace
    or Fourier transform.
    """
    # Boltzman constant
    kb = 1.38e-23
    # magnitude of the modulus
    G = (msd**-1.)*1./(np.pi*r*special.gamma(1.+alpha))
    G = kb*T*(1.e27)*G
    # Storage G1 and loss G2 moduli
    G1 = G*np.cos(np.pi*alpha/2.)
    G2 = G*np.sin(np.pi*alpha/2.)
    # Calculate omega
    omega = t**-1.
    # Convert omega to 1/s
    omega = omega*(1.e6)
    return [omega, G1, G2]


def shear_modulus_laplace_transform(t, msd, r, T, bw=0.01):
    """ Calculate the shear modulus by direct laplace transform of the MSD

    Parameters
    ----------
    t : 1d-array
        Time-lags in units of microseconds
    msd : 1d-array
          Mean-squared displacements at the time-lags
          ``t`` in units of nm^2
    r : float
        Radius of the probe particles in nanometers
    T : float
        Temperature in Kelvin
    bw : float, optional
         Bandwith parameter for analytic continuation of
         the laplace transform into fourier space by locally weighted
         regression

    Returns
    -------
    omega : 1d-array
            Vector of angular frequencies in units of 1/s
    G1 : 1d-array
         Storage modulus at angular frequencies ``omega`` in
         units of Pa
    G2 : 1d-array
         Loss modulus at angular frequencies ``omega`` in
         units of Pa

    Notes
    -----
    This function returns the shear modulus in the Fourier (omega) domain
    by first transforming to the Laplace (s) domain and then performing
    analytic continuation onto the imaginary axis in the complex plane.

    This method produces erroneous results in the upper and lower frequency
    extremes (roughly the first and last decade of frequencies).
    """
    # Boltzman constant
    kb = 1.38e-23
    # Calculate the direct laplace transform of the mean-squared displacement
    s = t**-1.
    msd_laplace = utils.laplace(t, msd, s)
    # Calculate the G*s in laplace space
    Gs = (msd_laplace**-1.)/(np.pi*r*s)
    # Perform a local power-law analysis of the laplace space shear modulus
    [Theta, logGs] = utils.loess(np.log(s), np.log(Gs), degree=1, alpha=bw)
    # Calculate the local scaling exponent
    alpha_direct = Theta[1, :]
    G = np.exp(logGs)
    G = kb*T*(1.e27)*G
    # Storage G1 and loss G2 moduli
    G1 = G*np.cos(np.pi*alpha_direct/2.)
    G2 = G*np.sin(np.pi*alpha_direct/2.)
    # Calculate omega
    omega = t**-1.
    # Convert omega to 1/s
    omega = omega*(1.e6)
    return [omega, G1, G2]


def full_dlsur_analysis(t, corr, ergodic, r, T, q, Ip, Ie,
                        calc_g1_kws={}, pwr_law_kws={}):
    """ Perform a full microrheology analysis from the correlation function.

    This function returns a table reporting particle motion statistics
    and the material shear modulus, beginning with the correlation
    coefficient. It implicitly calculates the intercept of the correlation
    function and a power-law analysis of the MSD to perform an algebraic
    calculation of the shear modulus. 

    Parameters
    ----------

    t : 1d-array
        Vector of lag-times (in microseconds)
    corr : 1d-array
           Vector of values of the correlation coefficient, corresponding
           to the time-lags ``t``.
    ergodic : boolean
              If ``ergodic==False`` then corrections are applied to the
              calculation of `g1`, based on ``Ip`` and ``Ie``
    r : float
        Particle radius in nanometers
    T : float
        Temperature in Kelvin
    q : float
        Scattering vector in 1/nm
    Ip : float, `optional`
         Scattering intensity at the meausurement position where ``corr``
         is collected
    Ie : 1d-array
        Vector of scattering intensities at different measurement positions
        in the cuvette
    calc_g1_kws : dictionary, `optional`
                  Dictionary of keyword arguments to pass to
                  ``analysis.tools.calc_g1()`` for computing the
                  intermediate scattering function
    pwr_law_kws : dictionary, `optional`
                  Dictionary of keyword arguments to pass to
                  ``analysis_tools.msd_local_pwr_law()`` for local
                  power-law analysis of the msd

    Returns
    -------
    dlsmicro_df : DataFrame
                  Dataframe containing table of results from DLS microrheology
                  analysis
    """

    # Find the intermediate scattering function
    g1 = calc_g1(t, corr, ergodic, Ip=Ip, Ie=Ie, **calc_g1_kws)

    # Calculate the power-law smoothing of the msd
    [msd_smooth, alpha] = msd_local_pwr_law(t, g1, q,
                                            **pwr_law_kws)

    # Calculate the shear modulus from the power-law smoothing
    [omega, G1, G2] = shear_modulus(t, msd_smooth, alpha, r, T)

    dlsmicro_dict = {'t': t, 'msd_smooth': msd_smooth,
                     'alpha': alpha, 'omega': omega,
                     'G1': G1, 'G2': G2}
    dlsmicro_df = pd.DataFrame(dlsmicro_dict)

    return dlsmicro_df
