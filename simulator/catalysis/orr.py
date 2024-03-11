# Catalysis Simulation Object
# made by Noh           2021. 8.
# modified by J. Park   2021.10. 2023.01
'''
solver catalysis module
run ORR, HER, and OER
'''
import os, sys
import numpy as np
import importlib
import json
from .catmodels import Catmodels
from .calgibbs import *
from .gibbsplot import *
from .aux import fix_atoms
from ...utils.units import T         # T is imported here from thermo.py
from ...utils.auxil import list2str, list2_format

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



### Workflow for the calculation of catalysis

def runORR(calc, sim_param , mode='opt', fix=None, act_site=None, vib=1, flabel='test', pH=0, cat_rxn='orr'):
    '''
    input
        T       not given from argument
                read from ...utils.units
        cat_rxn     ORR/OER
        act_site  atom index start from 0 to anchor adsorbate
    '''
    #print(f"0 act_site {act_site}")
    ### 1. DFT calculation for all the struct (DFT-opt & zpe)
    ### ORR and OER proceed DFT calculation in the same routine but returns different order
    totE, zpe, TS = run_series_ORR(calc, sim_param, mode, fix, act_site, vib, flabel, T, cat_rxn=cat_rxn)
    print(f"{'total energy':^15} {list2_format(totE)}\n{'zpe':^15} {list2_format(zpe)}\n{'Entropy':^15} {list2_format(TS)}")
    with open(f"{cat_rxn}_{flabel}.json", 'w') as f:
        f.write(json.dumps({'total energy': totE, 'zpe': zpe, 'TS': TS}))
    with open(f"{cat_rxn}_{flabel}.dat", 'w') as f:
        f.write( list2str(totE, decimal=2, head='totE', end='\n'))
        f.write( list2str( zpe, decimal=2, head='zpe',  end='\n'))
        f.write( list2str(  TS, decimal=2, head='TS',   end='\n'))
    ### 2. Gibbs energy calculation by reading OUTCAR: import analysis
    if vib:
        if cat_rxn == 'orr':
            Gibbs_pH   = calc_gibbs_ORR_4e_pH(totE=totE, zpe=zpe, TS=TS, pH=pH, Temp=T)
        elif cat_rxn == 'oer':
            Gibbs_pH   = calc_gibbs_OER_4e_pH(totE=totE, zpe=zpe, TS=TS, pH=pH, Temp=T)
    else:
        Gibbs_novib = calc_gibbs_ORR_4e_pH(totE=totE, pH=pH, Temp=T)
    #print(f"G_ORR: {Gibbs_novib}\nG_ORR_vib : {Gibbs_vib}")
    print(f"Gibbs_pH {list2_format(Gibbs_pH)}")
    
    ### 3. Plot Gibbs energy for the series of structures: import gibbsplot
    plot_OXR_4e_wU(Gibbs_pH, pH=pH, Temp=T, cat_rxn=cat_rxn)
    return totE, zpe, TS


def run_series_ORR(calc, sim_params, mode, fix, act_site, vib, flabel, Temp, cat_rxn):
    '''
    Make intermediates models and Do DFT calculation
    cat_rxn ORR and OER
    input
        mode        opt (default)
        fix         None (default)
                    b1L  fixed bottom 1 layer in case of slab
        flabel  
        act_site      passed to 

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
    fsuffix     = f"{flabel}_{irc}_cat"
    outfile     = f"{simulator.checkfile}_{fsuffix}"
    if not os.path.isfile(outfile):
        calc.run_calculator(mode=mode, fix=fixed_atoms)
        calc.save_files(fsuffix)
    totE_cat    = calc.get_total_energy(output_name=outfile)

    ### No vib cal for pure catalyst: vib for only adsorbate                                          

    ###### Make Intermediates geometry
    ### act_site is fixed here before atoms are disturbed
    #print(f"1: act_site {act_site}")
    if not act_site:
        act_site = calc.atoms.select_pivot(site='center')
    print(f"act_site index is {act_site}")
    interm_fnames   = ['O2', 'OOH', 'O', 'OH']
    
    catalyst_opt    = simmodule.readAtomicStructure(calc.optfile)
    orr_model       = Catmodels(catalyst_opt)

    interm_models   = orr_model.four_electron_intermediates_gen(mode='ORR', act_site=act_site)
        
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
        fsuffix = f"{flabel}_{i+1}_cat{interm_fnames[i]}"
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
