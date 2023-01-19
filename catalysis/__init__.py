#
# VASP Simulation Object
# made by Noh           2021. 8.
# modified by J. Park   2021.10. 2023.01

import os 
from .catmodels import Catmodeling
from .analysis import *
from .gibbsplot import *
from ..simulator.vasp import Vasp
from ..io import read_poscar

import operator
import sys
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
    Gibbs_novib = gibbs_HER([totE_sys], [totE_sysH])
    Gibbs_vib   = gibbs_HER([totE_sys], [totE_sysH], [zpe], [ts])

    ### 3. Plot Gibbs energy for the series of structures
    G_H_legend = ['noVib', 'Vib']
    G_H        = Gibbs_novib + Gibbs_vib

    plot_HER(G_H, G_H_legend)

    return 0

#def runORR(atoms, mode='opt', nproc=1, npar=1, encut=400, kpoints=[1,1,1],
#                   ediff = 0.0001, ediffg = -0.05, fix=None, active=None, vib=1, label='test'):
def runORR(calc, vas_param , nproc=1, mode='opt', fix=None, active=None, vib=1, label='test'):
    print(vas_param)
    ### 1. Run HER for a given system: generate several structures then calculate (opt & zpe)
    #tE, zpe, ts = run_series_ORR(atoms, mode=mode, nproc=nproc, npar=npar, encut=encut, kpoints=kpoints,\
    #                                ediff=ediff, ediffg=ediffg, fix=fix, active=active, vib=vib, label=label)
    tE, zpe, ts = run_series_ORR(calc, vas_param, nproc=nproc, mode=mode,fix=fix, active=active, vib=vib, label=label)
    print(f"total energy: {tE}\nZPE: {zpe}\nEntropy: {ts}")
    ### 2. Gibbs energy calculation by reading OUTCAR: import analysis
    #Gibbs_novib    = gibbs_ORR_4e_acid(TE=tE, pH=0)
    Gibbs_vib       = gibbs_ORR_4e_acid(TE=tE, ZPE=zpe, TS=ts)
    #print(f"G_ORR: {Gibbs_novib}\nG_ORR_vib : {Gibbs_vib}")
    print(f"G_ORR_vib : {Gibbs_vib}")
    ### --> find onset potential for input
    pot_onset       = - max(list(map(operator.sub, Gibbs_vib[1:],Gibbs_vib[:-1]))) 
    print(pot_onset)
    ### 3. Plot Gibbs energy for the series of structures: import gibbsplot
    plot_ORR_4e_acid(Gibbs_vib, U=pot_onset, legend=['U=1.23V', f'U={pot_onset:5.2f}V', 'U=0.00V'])
    return 0

def run_series_HER(atoms, mode='opt', nproc=1, npar=1, encut=400, kpoints=[1,1,1], 
                   ediff = 0.0001, ediffg = -0.05, fix=None, active=None, vib=1, label='test'):

   
    ### Model: substrate
    n_atoms = len(atoms)
    atoms_HER = Vasp(atoms)     ### 
    atoms_HER.run_VASP(mode=mode, nproc=nproc, npar=npar, encut=encut, kpoints=kpoints, \
                            ediff=ediff, ediffg=ediffg, fix=fix)
    
    os.system('mv OUTCAR OUTCAR_%s_sys' % label)
    os.system('mv XDATCAR XDATCAR_%s_sys' % label)

    TE_sys = atoms_HER.get_total_energy(output_name='OUTCAR_%s_sys' % label) ## out of class

    #from nanocore.io import read_poscar

    atoms_opt = read_poscar('CONTCAR')
    atoms2 = Catmodeling(atoms_opt) 
    atomsH = atoms2.HER_transition_gen(active=active)
    
    atomsH_HER = Vasp(atomsH)
    atomsH_HER.run_VASP(mode=mode, nproc=nproc, npar=npar, encut=encut, kpoints=kpoints, \
                             ediff=ediff, ediffg=ediffg, fix=fix)
    
    os.system('mv OUTCAR OUTCAR_%s_sysH' % label)
    os.system('mv XDATCAR XDATCAR_%s_sysH' % label)  

    TE_sysH = atomsH_HER.get_total_energy(output_name='OUTCAR_%s_sysH' % label)

    if vib:
        atomsH_opt = read_poscar('CONTCAR')
        
        fix_vib = []
        for i in range(n_atoms):
            idx = i+1
            fix_vib.append(idx)

        atomsH_Vib = Vasp(atomsH_opt)
        atomsH_Vib.run_VASP(mode='vib', nproc=nproc, npar=npar, encut=encut, kpoints=kpoints, \
                                 ediff=ediff, ediffg=ediffg, fix=fix_vib)

        os.system('mv OUTCAR OUTCAR_%s_sysH_Vib' % label)
        
        ZPE, TS = atomsH_Vib.get_vibration_energy(output_name='OUTCAR_%s_sysH_Vib' % label)

    if vib:
        return float(TE_sys), float(TE_sysH), float(ZPE), float(TS)
    else:
        return float(TE_sys), float(TE_sysH)

#def run_series_ORR(atoms, mode='opt', nproc=1, npar=1, encut=400, kpoints=[1,1,1],
#                ediff = 0.0001, ediffg = -0.05, fix=None, active=None, vib=1, label='test'):
def run_series_ORR(calc, vas_param, nproc=1, mode='opt', fix=None, active=None, vib=1, label='t'):
    '''
    Used parameter: mode, vib, label  
    Passed params : fix,  active
    '''
    print(f"in run_series_ORR: {vas_param}")
    irc = 0
    natoms          = len(calc.atoms)
    #orr_cat         = Vasp(atoms)
    ### Skip run_VASP if there is OUTCAR
    fsuffix         = f"{label}_{irc}_cat"
    outcar          = f"OUTCAR_{fsuffix}"
    print(f"outcar name: {outcar}")
    if not os.path.isfile(outcar):
        #orr_cat.run_VASP(mode=mode, nproc=nproc, npar=npar, encut=encut, kpoints=kpoints, \
        #                 ediff=ediff, ediffg=ediffg, fix=fix)
        calc.run_VASP(mode=mode, fix=fix)
        os.system(f'cp POSCAR  POSCAR_{fsuffix}')
        os.system(f'mv OUTCAR  {outcar}')
        os.system(f'mv XDATCAR XDATCAR_{fsuffix}')
        os.system(f'cp CONTCAR CONTCAR_{fsuffix}')
    totE_cat      = calc.get_total_energy(output_name=outcar)

    ### No vib cal for pure catalyst: vib for only adsorbate                                          

    ### Make Intermediates POSCAR
    INT_fnames      = ['O2', 'OOH', 'O', 'OH']  
    catalyst_opt    = read_poscar('CONTCAR')
    orr_model_gen   = Catmodeling(catalyst_opt)

    #orr_catO2, orr_catOOH, orr_catO, orr_catOH = orr_models.four_electron_transition_gen(mode='ORR', active=active)
    INT_models = orr_model_gen.four_electron_transition_gen(mode='ORR', active=active)
    
    #INT_models  = [orr_catO2, orr_catOOH, orr_catO, orr_catOH]
    
    ### list for data
    ltotE       = [totE_cat]
    lzpe        = [float(0.000)]
    lTS        = [float(0.000)]
    
    fix_vib      = []
    for j in range(natoms):
        idx = j+1
        fix_vib.append(idx)

    ### each model calculation
    print("start intermediates calculation"
    for i in range(len(INT_models)):
        calc = Vasp(INT_models[i])
        calc.set_options(**vas_param)

        ### skip if there is OUTCAR
        suffix = f"{label}_{i+1}_cat{INT_fnames[i]}"
        outcar = f"OUTCAR_{suffix}"
        if not os.path.isfile(outcar):
            #cal.run_VASP(mode=mode, nproc=nproc, npar=npar, encut=encut, kpoints=kpoints, \
            #     ediff=ediff, ediffg=ediffg, fix=fix)
            calc.run_VASP(mode=mode, fix=fix)
            os.system(f'cp POSCAR  POSCAR_{suffix}')
            os.system(f'mv OUTCAR  {outcar}')
            os.system(f'mv XDATCAR XDATCAR_{suffix}')  
            os.system(f'cp CONTCAR CONTCAR_{suffix}')  
        E = cal.get_total_energy(output_name=f'OUTCAR_{suffix}')
        ltotE.append(E)
        ### vib calculation for adsorbate
        if vib:
            suffix      += '_vib'
            outcar      = f"OUTCAR_{suffix}"
            opt_poscar  = read_poscar('CONTCAR')
            calc_vib    = Vasp(opt_poscar)
            calc_vib.set_options(**vas_param)
            if not os.path.isfile(outcar):
                #cal_vib.run_VASP(mode='vib', nproc=nproc, npar=npar, encut=encut, kpoints=kpoints, \
                #               ediff=ediff, ediffg=ediffg, fix=fix_vib)
                cal_vib.run_VASP(mode='vib', fix=fix_vib)
                os.system(f'mv OUTCAR {outcar}')
            zpe, TS = cal_vib.get_vibration_energy(output_name=f'{outcar}')
            lzpe.append(zpe)
            lTS.append(TS)
    if vib:
        return ltotE, lzpe, lTS
    else:
        return ltotE
