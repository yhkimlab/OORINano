# Catalysis Simulation Object
# made by Noh           2021. 8.
# modified by J. Park   2021.10. 2023.01
'''
solver catalysis module
run ORR, HER, and OER
'''
import os, sys
import numpy as np
import json
import importlib
from .catmodels import Catmodels
from .calgibbs import *
from .gibbsplot import *
from .aux import fix_atoms

from ...utils.units import T         # T is imported here from thermo.py

#...utils.auxililiary functions

### Workflow for the calculation of HER in catalysis

def runHER(calc, sim_params, atoms_list=None, mode='opt', fix=None, act_site=None, vib=1, flabel='test', pH=0):
    '''
    loop might be here for calc list
    '''
    simulator = calc.__class__
    simmodule = importlib.import_module(simulator.__module__)
    ### 1. Run HER for a given system: generate several structures then calculate (opt & zpe)
    if atoms_list:
        gibbsH=[]
        ### this doesnot make subdirectory for each models
        for i, atoms in enumerate(atoms_list):
            calc    = simulator(atoms)
            calc.set_options(**sim_params)
            print(f"calculate {flabel[i]}")
            totE_sys, totE_sysH, zpe, TS = run_series_HER(calc, sim_params, mode, fix, act_site, vib, flabel[i], T)
            gibbs_vib   = calc_gibbs_HER_pH(totE_sys, totE_sysH, zpe, TS, pH=pH, Temp=T)
            gibbsH.append(gibbs_vib)
            
        plot_HER(gibbsH, legend=flabel, pH=pH)

    else:
        totE_sys, totE_sysH, zpe, TS = run_series_HER(calc, sim_params, mode, fix, act_site, vib, flabel, T)
        print(f"total energy: {totE_sys} {totE_sysH}\nZPE: {zpe}\nEntropy: {TS}")

        
        ### 2. Gibbs energy calculation by reading OUTCAR
        Gibbs_novib = calc_gibbs_HER_pH(totE_sys, totE_sysH, pH=pH, Temp=T)
        Gibbs_vib   = calc_gibbs_HER_pH(totE_sys, totE_sysH, zpe, TS, pH=pH, Temp=T)
        print(f"G_HER: {Gibbs_novib}\nG_HER_vib : {Gibbs_vib}")
        ### 3. Plot Gibbs energy for the series of structures
        G_H_legend = ['Norskov', 'Vib']
        G_H        = [Gibbs_novib, Gibbs_vib]

        plot_HER(G_H, G_H_legend, pH=pH, Temp=T)
        #print(f"G_HER: {Gibbs_novib}\nG_HER_vib : {Gibbs_vib}")
    return 0

def run_series_HER(calc, sim_params, mode, fix, act_site, vib, flabel, Temp):
    '''
    mode        opt (default)
    fix         None (default)
                b1L  fixed bottom 1 layer in case of slab
    vib
    flabel  
    act_site      passed to 
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
    fsuffix     = f"{flabel}_cat"
    outfile     = f"{simulator.checkfile}_{fsuffix}"
    if not os.path.isfile(outfile):
        calc.run_calculator(mode=mode, fix=fixed_atoms)
        calc.save_files(fsuffix)
    totE_cat    = calc.get_total_energy(output_name=outfile)
    #print(f"totE_cat {totE_cat}")

    fsuffix     = f"{flabel}_catH"
    outfile     = f"{simulator.checkfile}_{fsuffix}"
   
    if not act_site:
        act_site = calc.atoms.select_pivot(site='center')
    #print(f"act_site {act_site}")

    catalyst_opt = simmodule.readAtomicStructure(calc.optfile)
    her_model = Catmodels(catalyst_opt) 
    atomsH = her_model.HER_intermediate_gen(act_site=act_site)
    
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
