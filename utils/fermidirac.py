'''
Common scientific functions such as Ferme-Dirac distribution
'''

import numpy as np
from scipy.misc import derivative
from ..utils.units import T, R

def fermi_dirac(E, fermi=0.0, t=T):
    return 1 / (1 + np.exp((E - fermi) / (t*R)))

def fd_derivative(E, t=T):
    fd_deriv = derivative(fermi_dirac, E, dx=1e-6)
    return fd_deriv

def fermi_dirac_farg(E, RT, fermi=0.0):
    '''
    
    RT      given from calling function of fd_derivative_weight
    '''
    return 1 / (1 + np.exp((E - fermi) / RT))

def fd_derivative_weight(E, weight_arg=0.4):
    '''
    This is weighting function of fermi_dirac_derivative with arbitarary parameter of RT
    weight_arg  required and passed to fermi_dirac_farg()
    '''

    fd_deriv = derivative(fermi_dirac_farg, E, dx=1e-6, args=(weight_arg,))
    return fd_deriv