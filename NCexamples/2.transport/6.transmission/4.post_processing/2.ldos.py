from nanocore import io, vis
from nanocore.simulator import siesta as s2
from matplotlib import cm
import numpy as np
import time, shutil,os,glob

cwd = os.getcwd()
os.chdir('../3.scatter+tbtrans/siesta')
files = os.listdir()
for f in files:
    if os.path.splitext(f)[-1] in ['.LDOS', '.XV']:
        label = os.path.splitext(f)[0]
        shutil.copy(f, cwd)
    elif f == 'STRUCT.fdf':
        shutil.copy(f, cwd)
os.chdir(cwd)

atoms = io.read_struct('STRUCT.fdf') 
cell = atoms.get_cell()
s2.get_ldos(cell[0], cell[1], cell[2], (0,0,0), (50, 50, 200), label)
os.system("xcrysden --xsf LDOS.XSF")