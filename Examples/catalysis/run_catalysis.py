import argparse
import re
import sys
from oorinano import catalysis
from oorinano.calculator.vasp import Vasp
from oorinano.calculator.vasp import readAtomicStructure as read_geo
from oorinano import surflab
from oorinano.utils.auxil import fname_ext, fname_root
import json

def run_catalysis(job, cat_rxn, pH, flabel, poscar, mode, Lvib, fix, nnode, nproc, Vasp_par ):
    '''
    job         run, model, incar, plot
    cat_rxn     ORR(default=orr), HER, OER
    flabel      filename of job to save OUTCAR, CONTCAR, XDATCAR (default='test')
    Loverwrite  Run or not(default) in case output file exists (OUTCAR_flabel_N_catO2[_vib] etc

    Run params:
        mode    opt sp vib
        Lvib    activate vibration calculation for adsorbate for TS (default=True)
        flabel  tag for filename to save OUTCAR, CONTCAR, XDATCAR, etc
        fix     None (default)
                b1L for slab to fix bottom 1 layer of slab
        pivot   adsorbate anchored position: list (position) or (int) atom index
                default: 1. the highest z, 2. center of xy-plane of supercell
    VASP params for INCAR, KPOINTS & mpirun:
        nproc
        kpoints
        poscar  Lvib
        INCAR:  mode Lvib Vasp_par keywords  npar, ispin, magmom, etc
            magmom will be list|dict such as ispin=2 & magmom=['N',2,...]|{'N':2,...(not checked)}
            mod change value in Vasp.wrtie_INCAR, Vasp.run_calculator
            ispin define here (below)
    '''
    
    ### 1. Make atoms
    
    if len(poscar) == 1:
        if re.search('pos', poscar[0], re.I) or re.search('con', poscar[0], re.I):
            ### read POSCAR
            atoms = read_geo(poscar[0])
        elif job == 'plot' and fname_ext(poscar[0]) == 'json':
            ### -p POSCAR is for input data file for plot
            fin = poscar[0]
    ## generate surface
    elif len(poscar) >= 3:
        if len(poscar) == 3:
            asize = int(poscar[2])
            size=(asize,asize,asize)
        elif len(poscar) == 5:
            size=(int(poscar[2]), int(poscar[3]), int(poscar[4]))
        print(f"fccsurf: {poscar[0]} {poscar[1]}, {size}")
        atoms = surflab.fccsurfaces(poscar[0], poscar[1], size, vac=15)
        ### for slab structure, fix bottom layer
        if fix is None:
            fix = 'b1L'
    else:
        atoms = None


    ### 2. Set params 
    ## npar vs ncore exclusive
    if not job == 'plot':
        if re.match('c', Vasp_par ):
            ncore=int(Vasp_par [1:])
        else:
            npar=int(Vasp_par [1:])

        if cat_rxn == 'orr' or cat_rxn == 'oer':
            
            ###     magmom = dict or list: ispin=2, magmom=['N',2]
            incar_params = dict(kpoints=[4,4,1], ediff=0.0001, ediffg=-0.05, encut=400, ispin=2)
            ### user incar test
            #incar_params['kband'] = "67 68 69 70"
            if 'npar' in locals():
                incar_params['npar']  = npar
            else:
                incar_params['ncore'] = ncore
            sim_params   = dict(nproc=nproc)
            sim_params.update(incar_params)
        elif cat_rxn == 'her':
            sim_params  = dict(npar=npar, kpoints=[1,1,1], nproc=nproc, ediff=0.01, ediffg=-0.04, encut=400, ispin=1)
        else:
            print(f"Error:: {cat_rxn} should be orr|oer|her")
            sys.exit(1)

        ### 3. Make Vasp instance and pass to runORR to make Vasp instance inside module 
        calc    = Vasp(atoms)
        calc.set_options(**sim_params)

    #calc.write_POSCAR(file_name='POSCAR.test')
    #sys.exit(10)

    ### 4. Run catalysis (VASP) | Show INCAR | Plot
    if job == 'run':
        if cat_rxn == 'orr' or cat_rxn ==  'oer':
            #
            catalysis.runORR(calc, sim_params, mode=mode, vib=Lvib, fix=fix, label=flabel, pH=pH, cat_rxn=cat_rxn)    #pivot = 24 (atom index)
        elif cat_rxn == 'her':
            catalysis.runHER(calc, sim_params, mode=mode, vib=Lvib, fix=fix, label=flabel)
    elif job == 'model':
        calc.write_POSCAR()
    elif job == 'incar':
        for k, v in calc.get_options():
            print(f"{k:>10}\t{v}")
    elif job == 'plot':
        ### read json file for data: totE, zpe, TS
        with open(fin) as f:
            ene = json.load(f)
        if re.search('orr', fin, re.I ):
            gibbs_vib = catalysis.calc_gibbs_ORR_4e(totE=ene['total energy'], ZPE=ene['zpe'], TS=ene['TS'], pH=pH, Temp=298.15)
            catalysis.plot_ORR_4e_wU(gibbs_vib, label=f'{fname_root(fin)}', pH=pH, Temp=298.15)
    else:
        print(f"Job error: select job in [run, model, incar, plot]")
        sys.exit(1)
    return 0

def main():
    parser = argparse.ArgumentParser(description="Running catalysis::\
                        \n\tselect catalytic job, subjob [run, show incar, ...], some options for vib, overwrite\
                        \n\tsystem params partition, node, etc are applied to specific system")
    parser.add_argument('-j', '--job', default='run', choices=['run', 'model','incar', 'plot'], help='incar: show default params')
    parser.add_argument('-c', '--cat_rxn', default='orr', choices=['orr', 'her', 'oer'], help='catalytic reactions')
    parser.add_argument('-ph', '--pH', default=0, type=int, help='get pH from command line')
    parser.add_argument('-l', '--flabel', default='test', help='label for dirname')
    parser.add_argument('-p', '--poscar', nargs='*', default=['POSCAR'], help="use any poscar or generate surface: ['Pt', '111', (3,3,3)]=[metal, surface index, size] or json file for job==plot")
    parser.add_argument('-t', '--test', action='store_true', help="change all the defaults")
    group_cat  = parser.add_argument_group(title='catalysis running options')
    group_cat.add_argument('-m', '--mode', default='opt', choices=['opt', 'sp'], help='Opt mode')
    group_cat.add_argument('-nv', '--novib', action='store_false', help='run orr without vibration')
    group_cat.add_argument('-fix', '--fix', help='fix bottom layer for slab')
    group_sys   = parser.add_argument_group(title='System-dependent inputs')
    group_sys.add_argument('-jn', '--jname', default='test', help='submit job name')
    group_sys.add_argument('-x', '--partition', help='partition name')
    group_sys.add_argument('-n', '--nnode', default=1, type=int, help='number of nodes: if needed')
    group_sys.add_argument('-np', '--nproc', type=int, default=24, help='number of process for mpirun')
    parallel = group_sys.add_mutually_exclusive_group()
    parallel.add_argument('--npar', type=int, default=4, help='npar value in INCAR')
    parallel.add_argument('--ncore', type=int, help='ncore value in INCAR for KISTI')
    parser.add_argument('-u', '--usage', action='store_true', help='explains how to run.')

    args = parser.parse_args()
    if args.usage:
        print(f"Usage::\
                \n    These are examples for job submit in queue systems and direct run\
                \n    Check 'readme.txt' to set VASP envirionment\
                \n    Run:\
                \n\t1. Submit jobscript with jobname, (slurm: partition, nnode, nproc) with variables\
                \n\t    slurm-cat::\
                \n\t\t$sbatch -J {args.jname} -p {args.partition} -N {args.nnode} -n {args.nproc}  slm_catalysis.sh\
                \n\t\t$sbatch -J {args.jname} -p {args.partition} -N {args.nnode} -n {args.nproc} --export=cat='{args.cat_rxn}',pos='gen' slm_catalysis.sh (default)\
                \n\t\t$sbatch -J {args.jname} -p {args.partition} -N {args.nnode} -n {args.nproc} --export=cat='{args.cat_rxn}',pos='cp' slm_catalysis.sh\
                \n\t\t    - pos='cp' for copy existing poscar & other char for slab generation\
                \n\t    pbs  ::\
                \n\t\t$qsub -N {args.jname} pbs_vasp_kisti_skl.sh\
                \n\t    Output::\
                \n\t\t/test     job directory is generated\
                \n\t\trun_catalysis.py is run inside job script\
                \n\t\tjob is running in subdir(jobname) & logfile is written in workdir (submit dir)\
                \n\t\tjob finishes: jobname.log -> jobname.out\
                \n\t2. Direct run inside job directory\
                \n\t    ORR with POSCAR\
                \n\t\t$python ../run_catalysis.py -j run -p Pt 111 3 -np {args.nproc} --npar $npar [--ncore $ncore]\
                \n\t3. Plot in the job directory\
                \n\t\t$python ../run_catalysis.py -c orr[oer] -ph 14\
            ")
        sys.exit(0)
    if args.ncore:
        nparallel='c'+str(args.ncore)
    else:
        nparallel='p'+str(args.npar)
    if args.test:
        args.mode = 'sp'
        args.novib = False

    run_catalysis(args.job, args.cat_rxn, args.pH, args.flabel, args.poscar, args.mode, args.novib, args.fix, args.nnode, args.nproc, nparallel)

if __name__ == "__main__":
    main()
