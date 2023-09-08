'''
executable environment: in KISTI
'''

import os
home = os.environ['HOME']
#
#  VASP
#
###     HOME
### VASP executable
vasp_calculator  = f"{home}/bin/skl_i18.0.3.2022.7_vasp_std" 

### supply POTCAR directory
vaspini = f"{home}/sandbox/pyvasp/vasppot/"

vasp_POTCAR_LDA  = vaspini + '1.POTPAW.LDA.54.RECOMMEND'
vasp_POTCAR_PBE  = vaspini + '2.POTPAW.PBE.54.RECOMMEND'
vasp_POTCAR_PW91 = vaspini + '2.POTPAW.PBE.54.RECOMMEND'        # connect to correct directory

#
# SIESTA
#

# Siesta calculator location
siesta_calculator = '[siesta-installed-location]/siesta-4.1-b3/Obj/siesta'

# Pseudopentiential files location
siesta_psf_location = ''

# Pyprojection files location
import os
temporary_location  = os.getcwd()
siesta_pyprojection = temporary_location + '/pyprojection/pyprojection.py'

# Siesta utilities location
siesta_util_location = '[siesta-installed-location]/siesta-4.1-b3/Util'
siesta_util_band = 'Bands/gnubands'
siesta_util_dos  = 'Eig2DOS/Eig2DOS'
siesta_util_pdos = 'Contrib/APostnikov/fmpdos'
siesta_util_rho  = 'Contrib/APostnikov/rho2xsf'
siesta_util_vh   = 'Macroave/Src/macroave'

