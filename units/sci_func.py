'''
Common scientific functions such as Ferme-Dirac distribution
'''

import numpy as np
from scipy.misc import derivative
from ..units import T, R

def fermi_dirac(E, fermi=0.0):
        return 1 / (1 + np.exp((E - fermi) / (T*R)))

def derivative_fd(E):
        deri_fd = derivative(fermi_dirac, E, dx=1e-6)
        return deri_fd