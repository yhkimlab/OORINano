#
# VASP Simulation Object
# made by Noh           2021. 8.
# modified by J. Park   2021.10.

import os 
from .catmodels import Catmodeling as Modeling
from .analysis import *
from .gibbsplot import *
from ..simulator.vasp import Vasp
from ..io import read_poscar

"""
attributes
    def runHER()
    def runORR()
    def run_series_HER()
    def run_series_ORR()
"""

### Workflow for the calculation of catalysis

def runHER(atoms, mode='opt', nproc=1, npar=1, encut=400, kpoints=[1,1,1],
                   ediff = 0.0001, ediffg = -0.05, fix=None, active=None, vib=1, label='test'):
    ### 1. Run HER for a given system: generate several structures then calculate (opt & zpe)
    totE_sys, totE_sysH, zpe, ts = run_series_HER(atoms, mode=mode, nproc=nproc, npar=npar, encut=encut, kpoints=kpoints,
                                    ediff=ediff, ediffg=ediffg, fix=fix, active=active, vib=vib, label=label)
    ### 2. Gibbs energy calculation by reading OUTCAR
    Gibbs_noVib = gibbs_HER([totE_sys], [totE_sysH])
    Gibbs_Vib   = gibbs_HER([totE_sys], [totE_sysH], [zpe], [ts])

    ### 3. Plot Gibbs energy for the series of structures
    G_H_legend = ['noVib', 'Vib']
    G_H        = Gibbs_noVib + Gibbs_Vib

    plot_HER(G_H, G_H_legend)

    return 0


def runORR(atoms, mode='opt', nproc=1, npar=1, encut=400, kpoints=[1,1,1],
                   ediff = 0.0001, ediffg = -0.05, fix=None, active=None, vib=1, label='test'):
    ### 1. Run HER for a given system: generate several structures then calculate (opt & zpe)
    tE, zpe, ts = run_series_ORR(atoms, mode=mode, nproc=nproc, npar=npar, encut=encut, kpoints=kpoints,
                                    ediff=ediff, ediffg=ediffg, fix=fix, active=active, vib=vib, label=label)
    print(f"total energy: {tE}\nZPE: {zpe}\nEntropy: {ts}")
    ### 2. Gibbs energy calculation by reading OUTCAR
    Gibbs_noVib = gibbs_ORR_4e_acid(TE=tE, pH=0)
    Gibbs_Vib   = gibbs_ORR_4e_acid(TE=tE, ZPE=zpe, TS=ts)
    print(f"G_ORR: {Gibbs_noVib}\nG_ORR_vib : {Gibbs_Vib}")
    ### 3. Plot Gibbs energy for the series of structures
    plot_ORR_4e_acid(Gibbs_Vib, U=0.7, legend=['U=1.23V', 'U=0.70V', 'U=0.00V'])

    return 0

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

    #from nanocore.io import read_poscar

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
    os.system('cp CONTCAR CONTCAR_%s_Sys' % label)
    TE_ORR_Sys   = ORR_Sys.get_total_energy(output_name='OUTCAR_%s_Sys' % label)
                                                                                             
    #from nanocore.io import read_poscar
                                                                                             
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

    ### each model calculation
    for i in range(len(cal_target)):
        cal = Vasp(cal_target[i])
        cal.run_VASP(mode=mode, nproc=nproc, npar=npar, encut=encut, kpoints=kpoints, \
                     ediff=ediff, ediffg=ediffg, fix=fix)
        os.system('mv OUTCAR OUTCAR_%s_Sys%s'   % (label, cal_name[i]))
        os.system('mv XDATCAR XDATCAR_%s_Sys%s' % (label, cal_name[i]))  
        os.system('cp CONTCAR CONTCAR_%s_Sys%s' % (label, cal_name[i]))  
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





