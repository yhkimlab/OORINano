# Catalysis Simulation Object
# made by Noh           2021. 8.
# modified by J. Park   2021.10. 2023.01
'''
solver catalysis module
run ORR, HER, and OER
'''
import os, sys
import numpy as np
from .catmodels import Catmodels
from .calgibbs import *
from .gibbsplot import *
import importlib
from ...utils.units import T         # T is imported here from thermo.py
from ...utils.auxil import list2_2f
import json
#...utils.auxililiary functions


OERinORR=['cat', 'OOH', 'O', 'OH']
OER_order=['cat', 'OH', 'O', 'OOH']

def oer_ordering(li):
    ### from ORR 4e calc, convert ordering to OER
    newl = []
    if len(li) == 4:
        for i in range(len(li)):
            newl.append(li[OERinORR.index(OER_order[i])])
    else:
        print("Algorithm error in oer_ordering()")
        sys.exit(111)
    return newl

def fix_atoms(atoms, fix):
    if fix == 'b1L':
        atoms.select_atoms("gid", 0)
    else:
        pass
    return None

### Workflow for the calculation of catalysis

def runHER(calc, sim_params, mode='opt', fix=None, ipivot=None, vib=1, label='test', pH=0):
    
    ### 1. Run HER for a given system: generate several structures then calculate (opt & zpe)
    totE_sys, totE_sysH, zpe, TS = run_series_HER(calc, sim_params, mode, fix, ipivot, vib, label, T)
    print(f"total energy: {totE_sys} {totE_sysH}\nZPE: {zpe}\nEntropy: {TS}")
    
    ### 2. Gibbs energy calculation by reading OUTCAR
    Gibbs_novib = calc_gibbs_HER([totE_sys], [totE_sysH], pH=pH, Temp=T)
    Gibbs_vib   = calc_gibbs_HER([totE_sys], [totE_sysH], [zpe], [TS], pH=pH, Temp=T)
    print(f"G_HER: {Gibbs_novib}\nG_HER_vib : {Gibbs_vib}")
    ### 3. Plot Gibbs energy for the series of structures
    G_H_legend = ['noVib', 'Vib']
    G_H        = Gibbs_novib + Gibbs_vib

    plot_HER(G_H, G_H_legend, pH=pH, Temp=T)

    return 0

def runORR(calc, sim_param , mode='opt', fix=None, ipivot=None, vib=1, label='test', pH=0, cat_rxn='orr'):
    '''
    input
        T       not given from argument
                read from ...utils.units
        cat_rxn     ORR/OER
        ipivot  atom index start from 0 to anchor adsorbate
    '''
    #print(f"0 ipivot {ipivot}")
    ### 1. DFT calculation for all the struct (DFT-opt & zpe)
    ### ORR and OER proceed DFT calculation in the same routine but returns different order
    totE, zpe, TS = run_series_ORR(calc, sim_param, mode, fix, ipivot, vib, label, T, cat_rxn=cat_rxn)
    print(f"{'total energy':^15} {list2_2f(totE)}\n{'ZPE':^15} {list2_2f(zpe)}\n{'Entropy':^15} {list2_2f(TS)}")
    with open(f"{cat_rxn}_{label}.json", 'w') as f:
        f.write(json.dumps({'total energy': totE, 'zpe': zpe, 'TS': TS}))
    with open(f"{cat_rxn}_{label}.dat", 'w') as f:
        #f.write(f"{'totE':5}{*totE,}\n{'zpe':5}{*zpe,}\n{'TS':5}{*TS,}")
        f.write(f"{'totE':5}{*totE,}\n{'zpe':5}{*zpe,}\n{'TS':5}{*TS,}")
    ### 2. Gibbs energy calculation by reading OUTCAR: import analysis
    if vib:
        if cat_rxn == 'orr':
            Gibbs_vib_pH   = calc_gibbs_ORR_4e_pH(totE=totE, ZPE=zpe, TS=TS, pH=pH, Temp=T)
        elif cat_rxn == 'oer':
            Gibbs_vib_pH   = calc_gibbs_OER_4e_pH(totE=totE, ZPE=zpe, TS=TS, pH=pH, Temp=T)
    else:
        Gibbs_novib = calc_gibbs_ORR_4e_pH(totE=totE, pH=pH, Temp=T)
    #print(f"G_ORR: {Gibbs_novib}\nG_ORR_vib : {Gibbs_vib}")
    print(f"Gibbs_vib_pH {list2_2f(Gibbs_vib_pH)}")
    
    ### 3. Plot Gibbs energy for the series of structures: import gibbsplot
    if cat_rxn == 'orr':
        plot_ORR_4e_wU(Gibbs_vib_pH, label=f'ORR_4e_pH{pH}', pH=pH, Temp=T)
    elif cat_rxn == 'oer':
        plot_OER_4e_wU(Gibbs_vib_pH,label=f'OER_4e_pH{pH}', pH=pH, Temp=T)
    return totE, zpe, TS


def run_series_HER(calc, sim_params, mode, fix, pivot, vib, label, Temp):
    '''
    mode        opt (default)
    fix         None (default)
                b1L  fixed bottom 1 layer in case of slab
    vib
    label  
    pivot      passed to 
    '''
    simulator = calc.__class__
    simmodule = importlib.import_module(simulator.__module__)
    ### Model: substrate
    natoms = len(calc.atoms)
    
    ### fixed start from 0 ?
    ngroup = calc.atoms.make_groups()
    if fix:
        fix_atoms(calc.atoms, fix)
        fixed_atoms = calc.atoms.get_selected()
        #print(f"{fix} atoms fixed: {fixed_atoms}")
    else:
        fixed_atoms = None
        #print("No atoms are fixed")
        
    ### skip if there is calc.checkfile
    fsuffix     = f"{label}_cat"
    outfile     = f"{simulator.checkfile}_{fsuffix}"
    if not os.path.isfile(outfile):
        calc.run_calculator(mode=mode, fix=fixed_atoms)
        calc.save_files(fsuffix)
    totE_cat    = calc.get_total_energy(output_name=outfile)
    #print(f"totE_cat {totE_cat}")

    fsuffix     = f"{label}_catH"
    outfile     = f"{simulator.checkfile}_{fsuffix}"
   
    if not pivot:
        pivot = calc.atoms.select_pivot(site='center')
    print(f"pivot {pivot}")

    catalyst_opt = simmodule.readAtomicStructure(calc.optfile)
    her_model = Catmodels(catalyst_opt) 
    atomsH = her_model.HER_intermediate_gen(pivot=pivot)
    
    calc = simulator(atomsH)
    calc.set_options(**sim_params)
    if not os.path.isfile(outfile):
        calc.run_calculator(mode=mode, fix=fixed_atoms)
        calc.save_files(fsuffix)
        
    totE_catH = calc.get_total_energy(output_name=outfile)

    if vib:
        fsuffix += '_vib'
        outfile  = f"{simulator.checkfile}_{fsuffix}"
        optfile  = simmodule.readAtomicStructure(calc.optfile)

        suffixv     = '_vib'
        outfile     += suffixv
     
        fix_vib = []
        for i in range(natoms):
            idx = i+1
            fix_vib.append(idx)

        calc = simulator(optfile)
        calc.set_options(**sim_params)
        if not os.path.isfile(outfile):
            calc.run_calculator(mode='vib', fix=fix_vib)
            calc.save_checkfile(fsuffix)
        
        zpe, TS = calc.get_vibration_energy(output_name=outfile, Temp=Temp)

        return float(totE_cat), float(totE_catH), float(zpe), float(TS)
    else:
        return float(totE_cat), float(totE_catH)

def run_series_ORR(calc, sim_params, mode, fix, ipivot, vib, label, Temp, cat_rxn):
    '''
    Make intermediates models and Do DFT calculation
    ORR and OER
    input
        mode        opt (default)
        fix         None (default)
                    b1L  fixed bottom 1 layer in case of slab
        vib
        label  
        ipivot      passed to 

    return: different order and length for ORR(5) and OER(4)
        ltotE, lzpe, lTS
    '''
    simulator = calc.__class__
    simmodule = importlib.import_module(simulator.__module__)
    #print(f"in run_series_ORR: {sim_params}")
    irc = 0
    natoms      = len(calc.atoms._atoms)

    ### fixed start from 0
    ngroup = calc.atoms.make_groups()
    if fix:
        fix_atoms(calc.atoms, fix)
        fixed_atoms = calc.atoms.get_selected()
        #print(f"{fix} atoms fixed: {fixed_atoms}")
    else:
        fixed_atoms = None
        #print("No atoms are fixed")
    
    ### Skip run_calculator if there is __outfile of the simulator
    fsuffix     = f"{label}_{irc}_cat"
    outfile     = f"{simulator.checkfile}_{fsuffix}"
    if not os.path.isfile(outfile):
        calc.run_calculator(mode=mode, fix=fixed_atoms)
        calc.save_files(fsuffix)
    totE_cat    = calc.get_total_energy(output_name=outfile)

    ### No vib cal for pure catalyst: vib for only adsorbate                                          

    ###### Make Intermediates geometry
    ### pivot is fixed here before atoms are disturbed
    #print(f"1: pivot {pivot}")
    if not ipivot:
        ipivot = calc.atoms.select_pivot(site='center')
    print(f"pivot index is {ipivot}")
    interm_fnames   = ['O2', 'OOH', 'O', 'OH']
    
    catalyst_opt    = simmodule.readAtomicStructure(calc.optfile)
    orr_model       = Catmodels(catalyst_opt)

    interm_models   = orr_model.four_electron_intermediates_gen(mode='ORR', pivot=ipivot)
        
    ### list for data
    ltotE       = [totE_cat]
    lzpe        = [float(0.000)]
    lTS         = [float(0.000)]
    
    ### all the atoms of cat (natoms) are fixed starting index from 1
    fix_vib      = []
    for j in range(natoms):
        idx = j+1
        fix_vib.append(idx)
    
    ### DFT calculation for each model 
    for i in range(len(interm_models)):
        ### if OER, skip calc of catOO but keep index for output files
        if i==0 and cat_rxn == 'oer':
            continue
        #print(f"Make a new class instance with {i+1}th intermediates")
        calc = simulator(interm_models[i])
        calc.set_options(**sim_params)

        ### skip if there is calc.checkfile
        print(f"{len(interm_models[i])}")
        fsuffix = f"{label}_{i+1}_cat{interm_fnames[i]}"
        outfile = f"{simulator.checkfile}_{fsuffix}"
        if not os.path.isfile(outfile):
            calc.run_calculator(mode=mode, fix=fixed_atoms)
            calc.save_files(fsuffix)
        E  = calc.get_total_energy(output_name=outfile)
        
        ltotE.append(E)
        ### vib calculation for adsorbate
        #print("start vib calculation")
        if vib:
            fsuffix += '_vib'
            outfile  = f"{simulator.checkfile}_{fsuffix}"
            optfile  = simmodule.readAtomicStructure(calc.optfile)
            calc = simulator(optfile)
            calc.set_options(**sim_params)
            if not os.path.isfile(outfile):
                calc.run_calculator(mode='vib', fix=fix_vib)
                calc.save_checkfile(fsuffix)
                
            zpe, TS = calc.get_vibration_energy(output_name=f'{outfile}', Temp=Temp)
            lzpe.append(zpe)
            lTS.append(TS)
    if vib:
        ### different length for ORR and OER but the same order
        if cat_rxn == 'orr':
            return ltotE, lzpe, lTS
        elif cat_rxn == 'oer':
            return oer_ordering(ltotE), oer_ordering(lzpe), oer_ordering(lTS)
    else:
        if cat_rxn == 'orr':
            return ltotE
        else:
            return oer_ordering(ltotE)
