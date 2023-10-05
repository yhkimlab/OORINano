'''
executable environment: Home cluster
setting variables:

    VASP                SIESTA
    vasp_calculator     siesta_calculator
    vasp_POTCAR_LDA     siesta_default_location
    vasp_POTCAR_PBE
'''    
###    HOME
vasp_odir="/TGM/Apps/VASP/OLD_BIN/5.4.4/O2/NORMAL"
vasp_ndir="/TGM/Apps/VASP/5.4.4.pl2"
vasp_dirv6="/TGM/Apps/VASP/bin/6.3.1"

vasp_dir = vasp_odir

### VASP executable 
vasp_calculator  = vasp_dir + '/vasp.5.4.4.pl2.O2.NORMAL.std.x' 

### supply POTCAR directory
vaspini = '/TGM/Apps/VASP/POTCAR/'

vasp_POTCAR_LDA  = vaspini + '1.POTPAW.LDA.54.RECOMMEND'
vasp_POTCAR_PBE  = vaspini + '2.POTPAW.PBE.54.RECOMMEND'
vasp_POTCAR_PW91 = vaspini + '2.POTPAW.PBE.54.RECOMMEND'        # connect to correct directory

### SIESTA e

siesta_dir = '/home2/littleyu/opt/siesta-v4.1-b4'
### Siesta calculator location
siesta_calculator = siesta_dir + '/Obj/siesta'

### Default siesta location
siesta_default_location='/home/joonho/pymod/nanocore/simulators/siesta/siesta_default'
### psf, fdf locates under default siesta locationunder with subdirectories of /psf, /model/elec, /model/channel
#siesta_psf_location = siesta_default_location + '/psf'


# Pyprojection files location
import os
temporary_location  = os.getcwd()
siesta_pyprojection = temporary_location + '/pyprojection/pyprojection.py'

# Siesta utilities location
siesta_util_tbtrans = siesta_dir + '/Util/TS/TBtrans/tbtrans'
siesta_util_band = 'Util/Bands/gnubands'
siesta_util_dos  = 'Util/Eig2DOS/Eig2DOS'
siesta_util_pdos = 'Util/Contrib/APostnikov/fmpdos'
siesta_util_rho  = 'Util/Contrib/APostnikov/rho2xsf'
siesta_util_vh   = 'Util/Macroave/Src/macroave'

