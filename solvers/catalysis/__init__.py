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
from ...io_1 import read_poscar

import operator
import sys

### Workflow for the calculation of catalysis

def runHER(calc, sim_params, mode='opt', fix=None, active=None, vib=1, label='test'):
    if calc.__module__.split('.')[-1] == 'vasp':
        read_geo = read_poscar
    else:
        print(f"read function error in simulator {calc.__module__}")
        sys.exit(11)
    ### 1. Run HER for a given system: generate several structures then calculate (opt & zpe)
    totE_sys, totE_sysH, zpe, TS = run_series_HER(calc, sim_params, read_geo, mode, fix, active, vib, label)
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
    ### define input function name depending on calc
    if calc.__module__.split('.')[-1] == 'vasp':
        read_geo = read_poscar
    else:
        print(f"read function error in simulator {calc.__module__}")
        sys.exit(11)
    print(sim_param)
    ### 1. Run ORR for a given system: generate several structures then calculate (opt & zpe)
    totE, zpe, TS = run_series_ORR(calc, sim_param, read_geo, mode, fix, active, vib, label)
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

def run_series_HER(calc, sim_params, read_geo, mode, fix, active, vib, label):
    '''
    read_geo    read function is required depending on simulator
    mode        opt (default)
    fix         None (default)
                b1L  fixed bottom 1 layer in case of slab
    vib
    label  
    active      passed to 
    '''
    simulator = calc.__class__
    ### Model: substrate
    natoms = len(calc.atoms)
    ### skip if there is calc.checkfile
    fsuffix     = f"{label}_cat"
    outfile     = f"{simulator.checkfile}_{fsuffix}"
    if not os.path.isfile(outfile):
        calc.run_catalysis(mode=mode, fix=fix)
        calc.save_files(fsuffix)
    totE_cat    = calc.get_total_energy(output_name=outfile)
    #print(f"totE_cat {totE_cat}")

    fsuffix     = f"{label}_catH"
    outfile     = f"{simulator.checkfile}_{fsuffix}"
   
    catalyst_opt = read_geo(calc.optfile)
    her_model = Catmodeling(catalyst_opt) 
    atomsH = her_model.HER_transition_gen(active=active)
    
    calc = simulator(atomsH)
    calc.set_options(**sim_params)
    if not os.path.isfile(outfile):
        calc.run_catalysis(mode=mode, fix=fix)
        calc.save_files(fsuffix)
        
    totE_catH = calc.get_total_energy(output_name=outfile)

    if vib:
        fsuffix += '_vib'
        outfile  = f"{simulator.checkfile}_{fsuffix}"
        optfile  = read_geo(calc.optfile)

        suffixv     = '_vib'
        outfile     += suffixv
        atomsH_opt  = read_geo(calc.optfile)
        
        fix_vib = []
        for i in range(natoms):
            idx = i+1
            fix_vib.append(idx)

        calc = simulator(optfile)
        calc.set_options(**sim_params)
        if not os.path.isfile(outfile):
            calc.run_catalysis(mode='vib', fix=fix_vib)
            calc.save_checkfile(fsuffix)
        
        zpe, TS = calc.get_vibration_energy(output_name=outfile)

        return float(totE_cat), float(totE_catH), float(zpe), float(TS)
    else:
        return float(totE_cat), float(totE_catH)

def run_series_ORR(calc, sim_params, read_geo, mode, fix, active, vib, label):
    '''
    read_geo    read function is required depending on simulator
    mode        opt (default)
    fix         None (default)
                b1L  fixed bottom 1 layer in case of slab
    vib
    label  
    active      passed to 
    '''
    simulator = calc.__class__
    #print(f"in run_series_ORR: {sim_params}")
    irc = 0
    natoms      = len(calc.atoms._atoms)

    ### fixed start from 0
    if fix == 'b1L':
        ngroup = calc.atoms.make_groups()
        calc.atoms.select_atoms("gid", 0)
        fixed = calc.atoms._selected
    elif type(fix) == list:
        fixed = fixed
    else:
        fixed = None
    ### Skip run_catalysis if there is __outfile of the simulator
    fsuffix     = f"{label}_{irc}_cat"
    outfile     = f"{simulator.checkfile}_{fsuffix}"
    if not os.path.isfile(outfile):
        calc.run_catalysis(mode=mode, fix=fixed)
        calc.save_files(fsuffix)
    totE_cat    = calc.get_total_energy(output_name=outfile)

    ### No vib cal for pure catalyst: vib for only adsorbate                                          

    ###### Make Intermediates geometry
    interm_fnames   = ['O2', 'OOH', 'O', 'OH']
    catalyst_opt    = read_geo(calc.optfile)
    orr_model       = Catmodeling(catalyst_opt)

    interm_models   = orr_model.four_electron_transition_gen(mode='ORR', active=active)
        
    ### list for data
    ltotE       = [totE_cat]
    lzpe        = [float(0.000)]
    lTS         = [float(0.000)]
    
    ### all the atoms of cat (natoms) are fixed starting index from 1
    fix_vib      = []
    for j in range(natoms):
        idx = j+1
        fix_vib.append(idx)
    
    ### each model calculation
    for i in range(len(interm_models)):
        #print(f"Make a new class instance with {i+1}th intermediates")
        calc = simulator(interm_models[i])
        calc.set_options(**sim_params)

        ### skip if there is calc.checkfile
        print(f"{len(interm_models[i])}")
        fsuffix = f"{label}_{i+1}_cat{interm_fnames[i]}"
        outfile = f"{simulator.checkfile}_{fsuffix}"
        if not os.path.isfile(outfile):
            calc.run_catalysis(mode=mode, fix=fixed)
            calc.save_files(fsuffix)
        E  = calc.get_total_energy(output_name=outfile)
        
        ltotE.append(E)
        ### vib calculation for adsorbate
        #print("start vib calculation")
        if vib:
            fsuffix += '_vib'
            outfile  = f"{simulator.checkfile}_{fsuffix}"
            optfile  = read_geo(calc.optfile)
            calc = simulator(optfile)
            calc.set_options(**sim_params)
            if not os.path.isfile(outfile):
                calc.run_catalysis(mode='vib', fix=fix_vib)
                calc.save_checkfile(fsuffix)
                
            zpe, TS = calc.get_vibration_energy(output_name=f'{outfile}')
            lzpe.append(zpe)
            lTS.append(TS)
    if vib:
        return ltotE, lzpe, lTS
    else:
        return ltotE
