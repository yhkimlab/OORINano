# Catalysis Simulation Object
# made by Noh           2021. 8.
# modified by J. Park   2021.10. 2023.01
'''
solver catalysis module
run ORR, HER, and OER
'''
import os 
from .catmodels import Catmodeling
from .analysis import *
from .gibbsplot import *
from ...io import read_poscar

import operator
import sys

### Workflow for the calculation of catalysis

def runHER(calc, sim_params, mode='opt', fix=None, active=None, vib=1, label='test'):
    ### 1. Run HER for a given system: generate several structures then calculate (opt & zpe)
    totE_sys, totE_sysH, zpe, TS = run_series_HER(calc, sim_params, mode, fix, active, vib, label)
    print(f"total energy: {totE_sys} {totE_sysH}\nZPE: {zpe}\nEntropy: {TS}")
    ### 2. Gibbs energy calculation by reading OUTCAR
    Gibbs_novib = gibbs_HER([totE_sys], [totE_sysH])
    Gibbs_vib   = gibbs_HER([totE_sys], [totE_sysH], [zpe], [TS])
    print(f"G_HER: {Gibbs_novib}\nG_HER_vib : {Gibbs_vib}")
    ### 3. Plot Gibbs energy for the series of structures
    G_H_legend = ['noVib', 'Vib']
    G_H        = Gibbs_novib + Gibbs_vib

    plot_HER(G_H, G_H_legend)

    return 0

def runORR(calc, sim_param , mode='opt', fix=None, active=None, vib=1, label='test'):
    print(sim_param)
    ### 1. Run ORR for a given system: generate several structures then calculate (opt & zpe)
    totE, zpe, TS = run_series_ORR(calc, sim_param, mode, fix, active, vib, label)
    print(f"total energy: {totE}\nZPE: {zpe}\nEntropy: {TS}")
    ### 2. Gibbs energy calculation by reading OUTCAR: import analysis
    Gibbs_novib    = gibbs_ORR_4e(totE=totE, pH=0)
    Gibbs_vib       = gibbs_ORR_4e(totE=totE, ZPE=zpe, TS=TS)
    print(f"G_ORR: {Gibbs_novib}\nG_ORR_vib : {Gibbs_vib}")
    ### --> find onset potential for input
    pot_onset       = - max(list(map(operator.sub, Gibbs_vib[1:],Gibbs_vib[:-1]))) 
    print(pot_onset)
    ### 3. Plot Gibbs energy for the series of structures: import gibbsplot
    plot_ORR_4e_acid(Gibbs_vib, U=pot_onset, legend=['U=1.23V', f'U={pot_onset:5.2f}V', 'U=0.00V'])
    return 0

def run_series_HER(calc, sim_params, mode, fix, active, vib, label):

    sim_class = calc.__class__
    ### Model: substrate
    natoms = len(calc.atoms)
    ### skip if there is OUTCAR
    fsuffix = f"{label}_cat"
    outcar = f"OUTCAR_{fsuffix}"
    #if not os.path.isfile(outcar):
    calc.run_catalysis(mode=mode, fix=fix)
    os.system(f'cp POSCAR  POSCAR_{fsuffix}')
    os.system(f'mv OUTCAR {outcar}')
    os.system(f'mv XDATCAR XDATCAR_{fsuffix}')
    os.system(f'cp CONTCAR CONTCAR_{fsuffix}')

    totE_cat = calc.get_total_energy(output_name=f'{outcar}') ## out of class
    #print(f"totE_cat {totE_cat}")

    fsuffix = f"{label}_catH"
    outcar  = f"OUTCAR_{fsuffix}"
    #if not os.path.isfile(outcar):
    catalyst_opt = read_poscar('CONTCAR')
    her_model = Catmodeling(catalyst_opt) 
    atomsH = her_model.HER_transition_gen(active=active)
    
    calc = sim_class(atomsH)
    calc.set_options(**sim_params)
    calc.run_catalysis(mode=mode, fix=fix)
    os.system(f'mv OUTCAR {outcar}')
    os.system(f'mv XDATCAR XDATCAR_{fsuffix}')  

    totE_catH = calc.get_total_energy(output_name=f'{outcar}')

    if vib:
        suffixv     = '_vib'
        outcar      += suffixv
        atomsH_opt  = read_poscar('CONTCAR')
        
        fix_vib = []
        for i in range(natoms):
            idx = i+1
            fix_vib.append(idx)

        calc = sim_class(atomsH_opt)
        calc.set_options(**sim_params)
        calc.run_catalysis(mode='vib', fix=fix_vib)

        os.system(f'mv OUTCAR {outcar}')
        
        zpe, TS = calc.get_vibration_energy(output_name=f'{outcar}')

    if vib:
        return float(totE_cat), float(totE_catH), float(zpe), float(TS)
    else:
        return float(TotE_cat), float(totE_catH)

def run_series_ORR(calc, sim_params, mode, fix, active, vib, label):
    '''
    mode    opt (default)
    fix     None (default)
            1L  fixed bottom 1 layer in case of slab
    vib
    label  
    active  passed to 
    '''
    sim_class = calc.__class__
    #print(f"in run_series_ORR: {sim_params}")
    irc = 0
    natoms      = len(calc.atoms)

    ### if fix=NL, sort atoms and fix NL

    ### Skip run_VASP if there is OUTCAR
    fsuffix     = f"{label}_{irc}_cat"
    outcar      = f"OUTCAR_{fsuffix}"
    if not os.path.isfile(outcar):
        calc.run_catalysis(mode=mode, fix=fix)
        os.system(f'cp POSCAR  POSCAR_{fsuffix}')
        os.system(f'mv OUTCAR  {outcar}')
        os.system(f'mv XDATCAR XDATCAR_{fsuffix}')
        os.system(f'cp CONTCAR CONTCAR_{fsuffix}')
        os.system(f'cp OSZICAR OSZICAR_{fsuffix}')
    totE_cat    = calc.get_total_energy(output_name=outcar)

    ### No vib cal for pure catalyst: vib for only adsorbate                                          

    ### Make Intermediates POSCAR
    interm_fnames   = ['O2', 'OOH', 'O', 'OH']  
    catalyst_opt    = read_poscar('CONTCAR')
    orr_model       = Catmodeling(catalyst_opt)

    interm_models = orr_model.four_electron_transition_gen(mode='ORR', active=active)
    #interm_models  = [orr_catO2, orr_catOOH, orr_catO, orr_catOH]
    
    ### list for data
    ltotE       = [totE_cat]
    lzpe        = [float(0.000)]
    lTS        = [float(0.000)]
    
    ### all the atoms of cat (natoms) are fixed starting index from 1
    fix_vib      = []
    for j in range(natoms):
        idx = j+1
        fix_vib.append(idx)
    
    ### each model calculation
    for i in range(len(interm_models)):
        #print(f"Make a new class instance with {i+1}th intermediates")
        calc = sim_class(interm_models[i])
        calc.set_options(**sim_params)

        ### skip if there is OUTCAR
        fsuffix = f"{label}_{i+1}_cat{interm_fnames[i]}"
        outcar = f"OUTCAR_{fsuffix}"
        if not os.path.isfile(outcar):
            calc.run_catalysis(mode=mode, fix=fix)
            os.system(f'cp POSCAR  POSCAR_{fsuffix}')
            os.system(f'mv OUTCAR  {outcar}')
            os.system(f'mv XDATCAR XDATCAR_{fsuffix}')  
            os.system(f'cp CONTCAR CONTCAR_{fsuffix}')
            os.system(f'cp OSZICAR OSZICAR_{fsuffix}')  # to check MAGMOM
        E = calc.get_total_energy(output_name=f'{outcar}')
        ltotE.append(E)
        ### vib calculation for adsorbate
        #print("start vib calculation")
        if vib:
            suffixv     = '_vib'
            outcar      += suffixv
            opt_poscar  = read_poscar('CONTCAR')
            calc = sim_class(opt_poscar)
            calc.set_options(**sim_params)
            if not os.path.isfile(outcar):
                calc.run_catalysis(mode='vib', fix=fix_vib)
                os.system(f'mv OUTCAR {outcar}')
            zpe, TS = calc.get_vibration_energy(output_name=f'{outcar}')
            lzpe.append(zpe)
            lTS.append(TS)
    if vib:
        return ltotE, lzpe, lTS
    else:
        return ltotE
