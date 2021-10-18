from ..atoms import *
from .. import io
from ..io import cleansymb, get_unique_symbs, convert_xyz2abc, ang2bohr
from ..units import ang2bohr
from glob import glob
import os, math
import numpy as np

from .catmodels import Catmodeling as Modeling
from .analysis import Analysis, Plotgibbs
from ..simulator.vasp import Vasp


#
# VASP Simulation Object
# made by Noh           2021. 8.
# modified by J. Park   2021.10.
#

def run_series_HER(atoms, mode='opt', nproc=1, npar=1, encut=400, kpoints=[1,1,1], 
                   ediff = 0.0001, ediffg = -0.05, fix=None, active=None, vib=1, label='test'):

   
    ### Model: substrate
    n_atoms = len(atoms)
    atoms_HER = Vasp(atoms)     ### 
    atoms_HER.run_VASP(mode=mode, nproc=nproc, npar=npar, encut=encut, kpoints=kpoints, \
                            ediff=ediff, ediffg=ediffg, fix=fix)
    
    os.system('mv OUTCAR OUTCAR_%s_Sys' % label)
    os.system('mv XDATCAR XDATCAR_%s_Sys' % label)

    TE_Sys = atoms_HER.get_total_energy(output_name='OUTCAR_%s_Sys' % label) ## out of class

    from NanoCore.io import read_poscar

    atoms_opt = read_poscar('CONTCAR')
    atoms2 = Modeling(atoms_opt) 
    atomsH = atoms2.HER_transition_gen(active=active)
    
    atomsH_HER = Vasp(atomsH)
    atomsH_HER.run_VASP(mode=mode, nproc=nproc, npar=npar, encut=encut, kpoints=kpoints, \
                             ediff=ediff, ediffg=ediffg, fix=fix)
    
    os.system('mv OUTCAR OUTCAR_%s_SysH' % label)
    os.system('mv XDATCAR XDATCAR_%s_SysH' % label)  

    TE_SysH = atomsH_HER.get_total_energy(output_name='OUTCAR_%s_SysH' % label)

    if vib:
        atomsH_opt = read_poscar('CONTCAR')
        
        fix_vib = []
        for i in range(n_atoms):
            idx = i+1
            fix_vib.append(idx)

        atomsH_Vib = Vasp(atomsH_opt)
        atomsH_Vib.run_VASP(mode='vib', nproc=nproc, npar=npar, encut=encut, kpoints=kpoints, \
                                 ediff=ediff, ediffg=ediffg, fix=fix_vib)

        os.system('mv OUTCAR OUTCAR_%s_SysH_Vib' % label)
        
        ZPE, TS = atomsH_Vib.get_vibration_energy(output_name='OUTCAR_%s_SysH_Vib' % label)

    if vib:
        return float(TE_Sys), float(TE_SysH), float(ZPE), float(TS)
    else:
        return float(TE_Sys), float(TE_SysH)

def run_series_ORR(atoms, mode='opt', nproc=1, npar=1, encut=400, kpoints=[1,1,1],                
                ediff = 0.0001, ediffg = -0.05, fix=None, active=None, vib=1, label='test'):
                                                                                             
 
    n_atoms      = len(atoms)
    ORR_Sys      = Vasp(atoms)
    ORR_Sys.run_VASP(mode=mode, nproc=nproc, npar=npar, encut=encut, kpoints=kpoints, \
                         ediff=ediff, ediffg=ediffg, fix=fix)
 
    os.system('mv OUTCAR OUTCAR_%s_Sys' % label)
    os.system('mv XDATCAR XDATCAR_%s_Sys' % label)
    TE_ORR_Sys   = ORR_Sys.get_total_energy(output_name='OUTCAR_%s_Sys' % label)
                                                                                             
    from NanoCore.io import read_poscar
                                                                                             
    System_opt   = read_poscar('CONTCAR')
    ORR_Sys_opt  = Modeling(System_opt)
    
    ORR_SysO2, ORR_SysOOH, ORR_SysO, ORR_SysOH = ORR_Sys_opt.four_electron_transition_gen(mode='ORR', active=active)
    
    #####

    cal_target   = [ORR_SysO2, ORR_SysOOH, ORR_SysO, ORR_SysOH]
    cal_name     = ['O2', 'OOH', 'O', 'OH']  
    
    TE           = [TE_ORR_Sys]
    E_ZPE        = [float(0.000)]
    E_TS         = [float(0.000)]
    
    fix_vib      = []
    for j in range(n_atoms):
        idx = j+1
        fix_vib.append(idx)

    for i in range(len(cal_target)):
        cal = Vasp(cal_target[i])
        cal.run_VASP(mode=mode, nproc=nproc, npar=npar, encut=encut, kpoints=kpoints, \
                     ediff=ediff, ediffg=ediffg, fix=fix)
        os.system('mv OUTCAR OUTCAR_%s_Sys%s'   % (label, cal_name[i]))
        os.system('mv XDATCAR XDATCAR_%s_Sys%s' % (label, cal_name[i]))  
        E = cal.get_total_energy(output_name='OUTCAR_%s_Sys%s' % (label, cal_name[i]))
        TE.append(E)

        if vib:
            cal_opt = read_poscar('CONTCAR')
            cal_vib = Vasp(cal_opt)
            cal_vib.run_VASP(mode='vib', nproc=nproc, npar=npar, encut=encut, kpoints=kpoints, \
                               ediff=ediff, ediffg=ediffg, fix=fix_vib)
            os.system('mv OUTCAR OUTCAR_%s_Sys%s_Vib' % (label, cal_name[i]))
            ZPE, TS = cal_vib.get_vibration_energy(output_name='OUTCAR_%s_Sys%s_Vib' % (label, cal_name[i]))
            E_ZPE.append(ZPE)
            E_TS.append(TS)
                                                                                             
    if vib:
        return TE, E_ZPE, E_TS
    else:
        return TE





