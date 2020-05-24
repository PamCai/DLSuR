import numpy as np


# stretched exponential decay function
def stretched_exp(x, x0, a, beta):
    return x0*np.exp(-a*(x**beta))


# 'exp exp' decay function
def expexp(x, a0, a1, lam, beta):
    return a0*np.exp(-a1*(1-np.exp(-lam*x**beta)))
