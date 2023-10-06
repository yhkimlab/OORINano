import argparse
import os
import re
import sys
from nanocore import catalysis
from nanocore.simulators.vasp import Vasp
from nanocore.simulators.vasp import readAtomicStructure as read_geo
from nanocore import surflab

def run_catalysis(job, cat_kind, flabel, Loverwrite, poscar, mode, Lvib, nnode, nproc, sparallel):
    '''
    job         [ORR(default=orr), HER, OER]
    subjob      [VASP run, show INCAR, plot Gibbs]
    flabel      filename of job to save OUTCAR, CONTCAR, XDATCAR (default='test')
    Loverwrite  Run or not(default) in case output file exists (OUTCAR_flabel_N_catO2[_vib] etc

    Run params:
        mode    ['opt'(default), 'sp'] for optimization or single point calculation
        vib     activate vibration calculation for adsorbate for TS (default=True)
        label   == flable, filename of job to save OUTCAR, CONTCAR, XDATCAR, etc
        fix     None (default)
                b1L for slab to fix bottom 1 layer of slab
        pivot   adsorbate anchored position: list (position) or (int) atom index
                default: 1. the highest z, 2. center of xy-plane of supercell
    VASP params for INCAR, KPOINTS & mpirun:
        nproc
        kpoints
        incar keywords  npar, ispin, magmom, etc
            magmom will be list|dict such as ispin=2 & magmom=['N',2,...]|{'N':2,...(not checked)} 
    '''
    
    ### 1. Make atoms
    ## read poscar
    if len(poscar) == 1 and re.search('pos', poscar[0], re.I):
        atoms = read_geo(poscar[0])
    ## generate surface
    elif len(poscar) >= 3:
        if len(poscar) == 3:
            asize = int(poscar[2])
            size=(asize,asize,asize)
        elif len(poscar) == 5:
            size=(int(poscar[2]), int(poscar[3]), int(poscar[4]))
        print(f"fccsurf: {poscar[0]} {poscar[1]}, {size}")
        atoms = surflab.fccsurfaces(poscar[0], poscar[1], size, vac=15)
    else:
        atoms = None


    ### 2. Set params 
    ## npar vs ncore exclusive
    if re.match('c', sparallel):
        ncore=int(sparallel[1:])
    else:
        npar=int(sparallel[1:])

    if cat_kind == 'orr':
        ### INCAR params: 
        ###     magmom = dict or list: ispin=2, magmom=['N',2]
        if 'npar' in locals():
            incar_params = dict(npar=npar, kpoints=[4,4,1], ediff=0.0001, ediffg=-0.05, encut=400, ispin=2)
        else:
            incar_params = dict(ncore=ncore, kpoints=[4,4,1], ediff=0.0001, ediffg=-0.05, encut=400, ispin=2)

        sim_params   = dict(nproc=nproc)
        sim_params.update(incar_params)
    elif cat_kind == 'her':
        sim_params  = dict(npar=npar, kpoints=[1,1,1], nproc=nproc, ediff=0.01, ediffg=-0.04, encut=400, ispin=1)
    else:
        print(f"Error:: {cat_kind} should be 'orr'|'her'")
        sys.exit(1)

    ### 3. Make Vasp instance and pass to runORR to make Vasp instance inside module 
    calc    = Vasp(atoms)
    calc.set_options(**sim_params)

    #calc.write_POSCAR(file_name='POSCAR.test')
    #sys.exit(10)

    ### 4. Run VASP | Show INCAR | Plot
    if job == 'run':
        if cat_kind == 'orr':
            catalysis.runORR(calc, sim_params, mode='opt', vib=True, fix='b1L', label=flabel, pH=14)    #pivot = 24 (atom index)
        elif cat_kind == 'her':
            catalysis.runHER(calc, sim_params, mode='opt', vib=False, fix='b1L', label=flabel)
    elif job == 'model':
        calc.write_POSCAR()
    elif job == 'incar':
        for k, v in calc.get_options():
            print(f"{k:>10}\t{v}")
    elif job == 'plot':
        ## To plot: get values of totE, zpe, TS and plot
        totE, zpe, TS = catalysis.runORR(calc, sim_params, mode=mode, vib=True, fix='b1L', label=flabel, pH=14)
        print(totE)
        print(zpe)
        print(TS)
        #catalysis.plot_ORR_4e_acid(G_ORR_vib, U=0.7, legend=['U=1.23V', 'U=0.70V', 'U=0.00V'])
    elif job == 'dos':
        pass
    else:
        pass
    return 0

def main():
    parser = argparse.ArgumentParser(description="Running catalysis::\
                        \n\tselect catalytic job, subjob [run, show incar, ...], some options for vib, overwrite\
                        \n\tsystem params partition, node, etc are applied to specific system")
    parser.add_argument('-j', '--job', default='run', choices=['run', 'model','incar', 'plot','test'], help='incar: show default params')
    parser.add_argument('-c', '--cat_kind', default='orr', choices=['orr', 'her', 'oer'], help='catalytic reactions')
    parser.add_argument('-l', '--flabel', default='test', help='label for dirname')
    parser.add_argument('-o', '--overwrite', action='store_true', help='if there exists dir, overwrite')
    parser.add_argument('-p', '--poscar', nargs='*', help="use any poscar or generate surface: ['Pt', '111', (3,3,3)]=[metal, surface index, size]")
    parser.add_argument('-t', '--test', action='store_true', help="change all the defaults")
    group_vasp  = parser.add_argument_group(title='NC running options')
    group_vasp.add_argument('-m', '--mode', default='opt', choices=['opt', 'sp'], help='Opt mode')
    group_vasp.add_argument('-nv', '--novib', action='store_false', help='run orr without vibration')
    group_sys   = parser.add_argument_group(title='System-dependent inputs')
    group_sys.add_argument('-N', '--nnode', default=1, type=int, help='number of nodes: if needed')
    group_sys.add_argument('-np', '--nproc', type=int, default=24, help='number of process for mpirun')
    parallel = group_sys.add_mutually_exclusive_group()
    parallel.add_argument('--npar', type=int, default=4, help='npar value in INCAR')
    parallel.add_argument('--ncore', type=int, help='ncore value in INCAR for KISTI')
    parser.add_argument('-u', '--usage', action='store_true', help='explains how to run.')

    args = parser.parse_args()
    if args.usage:
        print(f"Usage::\
                \n\tThis is 3 examples for job submit in queue systems and direct run\
                \n\tCheck 'readme.txt' to set VASP envirionment\
                \n\t    1. Run jobscript with jobname, (slurm: partition, nnode, nproc) with variables\
                \n\t\t(slurm): $sbatch -J test -p X3 -N 1 -n 20 --export=cat='{args.cat_kind}' slurm_sbatch_nc.sh\
                \n\t\t( pbs ): $qsub -N test pbs_vasp_kisti_skl.sh\
                \n\t\t    /test     job directory is generated\
                \n\t\t    run_catalysis.py is run inside job script\
                \n\t    2. Direct run inside job directory\
                \n\t\trun_catalysis.py -c orr -j run -np {args.nproc} [--npar $npar|--ncore $ncore]\
                \n\t** job is running in work dir(jobname) & logfile is written in submit dir\
                \n\t    mpirun runs in class Vasp\
                \n\t** job finishes: jobname.log -> jobname.out\
            ")
        sys.exit(0)
    if args.ncore:
        nparallel='c'+str(args.ncore)
    else:
        nparallel='p'+str(args.npar)
    if args.test:
        args.mode = 'sp'

    ### To test algorithm
    if args.job == 'test':
        args.job    = 'run'
        args.mode   = 'sp'
        args.novib  = True
        ### control more: ispin, kpoints

    run_catalysis(args.job, args.cat_kind, args.flabel, args.overwrite, args.poscar, args.mode, args.novib, args.nnode, args.nproc, nparallel)

if __name__ == "__main__":
    main()
