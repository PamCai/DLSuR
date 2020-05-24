import numpy as np
from numpy import linalg
from numpy import random
from scipy import integrate
from scipy.optimize import curve_fit

##########################################################
#gauusian_weight
#compute a gaussian weight matrix W, given a vector x
#the reference observation point i (index) and a bandwith
#tau
##########################################################

def gaussian_weight(x, i, tau):
    #compute the vector of weights
    w = np.exp(-((x-x[i])**2.)/(tau))
    #transform into a diagonal array
    W = np.diag(w)
    return W

#########################################################
#Weighted linaer regression
#perform linear regression given a matrix of inputs
#X, whose rows correspond to observations and whose columns
#correspond to features and a vector of outputs
#y
#inputs:
#X - matrix of inputs (mXn) with m features and n observations
#y - vector of ouput observations (size n)
#W - n X n weight matrix
#ouputs:
#theta - vector of best-fit parameters
#########################################################

def weighted_linear_reg(X,y,W):
    m1 = linalg.inv(np.dot(np.dot(X.T,W),X))
    m2 = np.dot(np.dot(X.T,W),y)

    theta = np.dot(m1,m2)
    return theta

#####################################################################
#Loess
#Compute a locally weighted linear regression given a vector
#of inputs x and a vector of outputs y
#Inputs
#x - vector of inputs
#y - vector of ouputs
#degree - degree of polynomial for fit (1 or 2)
#alpha - smoothing parameter. Defined as the ratio of tau**2. divided
#by the square in the difference between the two extrema in x


def loess(x,y,degree,alpha):
    n = len(x)
    #######################################################
    #First degree polynomial regression
    #######################################################
    if degree == 1:

        #First feature in the X array is a constant
        #Second feature is the value of x
        const = np.ones(n)
        X = np.array([const,x]).T
        W = np.identity(n)
        tau = alpha * np.sqrt((x[0]-x[-1])**2.)

        #Construct the theta matrix
        Theta = np.zeros([2,n])
        #Construct the prediction vector
        yp = np.zeros(n)
        #Loop over observations
        for i in range(n):
            #compute the weight matrix
            W = gaussian_weight(x,i,tau)
           #perform linear regression
            theta = weighted_linear_reg(X,y,W)
            Theta[:,i] = theta
            yp[i] = np.dot(X[i,:],theta)
        return [Theta,yp]

    #######################################################
    #Second degree polynomial regression
    #######################################################
    if degree == 2:

        #First feature in the X array is a constant
        #Second feature is the value of x
        #Third feature is the value of x**2.
        const = np.ones(n)
        X = np.array([const,x,x**2.]).T
        W = np.identity(n)
        tau = alpha * np.sqrt((x[0]-x[-1])**2.)

        #Construct the theta matrix
        Theta = np.zeros([3,n])

        #Construct the prediction vector
        yp = np.zeros(n)
        #Loop over observations
        for i in range(n):
            #compute the weight matrix
            W = gaussian_weight(x,i,tau)
            #perform linear regression
            theta = weighted_linear_reg(X,y,W)
            Theta[:,i] = theta
            yp[i] = np.dot(X[i,:],theta)
        return [Theta,yp]

##################################
#Numerical Laplace transform
##################################
#Inputs
#t - time vector
#f - vector of discrete values of function to be transformed
#S - vector of Laplace frequencies to evaluate the transform at


def laplace(t, f, S):
    L = np.zeros(len(S))
    for i, s in enumerate(S):
        y = f*np.exp(-s*t)
        L[i] = integrate.trapz(y, t)
    return L


def get_cross_validation_score(t, y, func, p0=None):
    """ Obtain the leave-one-out cross-validation score for
    a functional model of a data set over a given fitting window.

    Parameters
    ----------
    t : 1-d array
        Length N vector  of values for the independent variable of a dataset.
    y : 1-d array
        Length N vector of values for the dependent variable of a dataset.
    func : callable function
           Model to fit the data against. Must be of the form
           `f(t, p1, p2, ..., pM)` where `t` is the independent variable and
           `p1, p2, ..., pM` is a set of M parameters to fit.

    p0 : 1-d array or list, `optional`
         Initial guesses for the M parameters to fit, [p1, p2, ..., pM]

    Returns
    -------
    cv : float 
         Leave-one-out cross-validation score for the model
    """
    cv = 0.
    yskip = y
    tskip = t
    n = len(tskip)
    # Loop over each (t, y) point and perform a
    # fit to func with that point removed
    # calculate MSE for that point
    for i in range(n):
        ytest = np.delete(yskip, i)
        ttest = np.delete(tskip, i)
        try:
            paramsi = curve_fit(func, ttest, ytest, p0=p0, maxfev=10000)[0]
            yfiti = func(tskip[i], *paramsi)
            erri = (yskip[i]-yfiti)**2.
        except RuntimeError:
            erri = 1.e3
        cv = cv + erri
    # Average cv scores
    cv = cv/np.float(n)

    # If fit is not found for this window, penalize strongly
    try:
        paramsi = curve_fit(func, t, y, p0=p0, maxfev=10000)[0]
    except RuntimeError:
        cv = 1.e6
    return cv


def minimize_cv_error(t, y, twindows, func, p0=None):
    """ Find the fitting interval that minimizes the cross-validation
    error for a model fitted to a sub-interval of a dataset, 
    given a set of possible intervals in the independent variable 
    over which to perform the fit.

    Parameters
    ----------
    t : 1-d array
        Length N vector  of values for the independent variable of a dataset.
    y : 1-d array
        Length N vector of values for the dependent variable of a dataset.
    twindows : List of len 2 lists
               List of the form [[t0, tend1], [t0, tend2], ...]
               where [t0, tendi] represents the `ith` closed interval over which to
               fit the data set and find the cross-validation score.
    func : callable function
           Model to fit the data against. Must be of the form
           `f(t, p1, p2, ..., pM)` where `t` is the independent variable and
           `p1, p2, ..., pM` is a set of M parameters to fit.

    p0 : 1-d array or list, `optional`
         Initial guesses for the M parameters to fit, [p1, p2, ..., pM]

    Returns
    -------
    twindow_min : List
                  List of the form `[t0, tend]' where t0 and tend are the
                  beginning and end of the closed interval in ``t`` that
                  provides the lowest cross-validation score for the fitted
                  model
    pmin : list
           Optimal parameters for fitting the model ``func`` to the data
           over the sub-interval ``twindow_min``
    CV_min : float
             Leave-one-out cross-validation error for the model ``func``
             over the interval ``twindow_min``
    """
    CVs = []
    params = []
    for twindow in twindows:
        tinds = [np.argmin(np.abs(t-twindow[0])),
                 np.argmin(np.abs(t-twindow[1]))]
        tfit = t[tinds[0]: tinds[1]+1]
        yfit = y[tinds[0]: tinds[1]+1]
        CVs.append(get_cross_validation_score(tfit, yfit, func, p0))
        try:
            paramsi = curve_fit(func, tfit, yfit, p0=p0, maxfev=100000)[0]
            params.append(paramsi)
        except RuntimeError:
            params.append(None)
    CV_argmin = np.argmin(CVs)
    CV_min = CVs[CV_argmin]
    twindow_min = twindows[CV_argmin]
    pmin = params[CV_argmin]
    return [twindow_min, pmin, CV_min]


def laplace_merge(omega, G1, G2, G1_Fdirect, G2_Fdirect):
    """ Find the interval where the Laplace transform is valid over
    the frequency space and replace the Fourier transform data with 
    the Laplace transform modulus.

    Parameters
    ----------
    omega : 1-d array
            Vector of angular frequencies in units of rad/s
    G1 : 1-d array
         Storage modulus at angular frequencies ``omega`` in units of Pa
         calculated using Fourier transform
    G2 : 1-d array
         Loss modulus at angular frequencies ``omega`` in units of Pa
         calculated using Fourier transform
    G1_Fdirect: 1-d array
                Storage modulus at angular frequencies ``omega`` in units 
                of Pa calculated using direct Laplace transform
    G2_Fdirect: 1-d array
                Loss modulus at angular frequencies ``omega`` in units
                of Pa calculated using direct Laplace transform

    Returns
    -------
    G1_plot : 1-d array
              Loss modulus at angular frequencies ``omega`` in units of Pa
              after merging Laplace and Fourier transform moduli
    G2_plot : 1-d array
              Storage  modulus at angular frequencies ``omega`` in units of
              Pa after merging Laplace and Fourier transform moduli
    """
    G1_lower = 0
    G2_lower = 0
    G1_upper = len(omega) - 1
    G2_upper = len(omega) - 1
    for n in range(len(omega)-1):
        if (G1[n]-G1_Fdirect[n])/np.abs(G1[n]-G1_Fdirect[n]) != (G1[n+1]-G1_Fdirect[n+1])/np.abs(G1[n+1]-G1_Fdirect[n+1]):
            G1_lower = n
            break
    for n in range(len(omega)-1):
        if (G2[n]-G2_Fdirect[n])/np.abs(G2[n]-G2_Fdirect[n]) != (G2[n+1]-G2_Fdirect[n+1])/np.abs(G2[n+1]-G2_Fdirect[n+1]):
            G2_lower = n
            break
    for n in range(len(omega)-1):
        if (G1[len(omega)-1-n]-G1_Fdirect[len(omega)-1-n])/np.abs(G1[len(omega)-1-n]-G1_Fdirect[len(omega)-1-n]) != (G1[len(omega)-2-n]-G1_Fdirect[len(omega)-2-n])/np.abs(G1[len(omega)-2-n]-G1_Fdirect[len(omega)-2-n]):
            G1_upper = len(omega)-1-n
            break
    for n in range(len(omega)-1):
        if (G2[len(omega)-1-n]-G2_Fdirect[len(omega)-1-n])/np.abs(G2[len(omega)-1-n]-G2_Fdirect[len(omega)-1-n]) != (G2[len(omega)-2-n]-G2_Fdirect[len(omega)-2-n])/np.abs(G2[len(omega)-2-n]-G2_Fdirect[len(omega)-2-n]):
            G2_upper = len(omega)-1-n
            break

    G1_plot = np.zeros_like(G1)
    G2_plot = np.zeros_like(G2)
    G1_plot[0:G1_lower] = G1[0:G1_lower]
    G1_plot[G1_lower:G1_upper] = G1_Fdirect[G1_lower:G1_upper]
    G1_plot[G1_upper:] = G1[G1_upper:]
    G2_plot[0:G2_lower] = G2[0:G2_lower]
    G2_plot[G2_lower:G2_upper] = G2_Fdirect[G2_lower:G2_upper]
    G2_plot[G2_upper:] = G2[G2_upper:]
    return G1_plot, G2_plot


