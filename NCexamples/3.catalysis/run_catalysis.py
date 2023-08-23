import argparse
import os
import sys
from nanocore import io
from nanocore import catalysis
from nanocore.simulator.vasp import Vasp
from nanocore import surflab

def run_catalysis(job, cat_kind, flabel, Loverwrite, poscar, mode, Lvib, nnode, nproc, npar):
    '''
    job         [ORR(default=orr), HER, OER]
    subjob      [VASP run, show INCAR, plot Gibbs]
    flabel      filename of job to save OUTCAR, CONTCAR, XDATCAR (default='test')
    Loverwrite  Run or not(default) in case output file exists (OUTCAR_flabel_N_catO2[_vib] etc

    Run params:
        mode    ['opt'(default), 'sp'] for optimization or single point calculation
        vib     activate vibration calculation for adsorbate for TS (default=True)
        label   filename of job to save OUTCAR, CONTCAR, XDATCAR, etc
        fix     None (default)
                b1L for slab to fix bottom 1 layer of slab
        active
    VASP params for INCAR, KPOINTS & mpirun:
        nproc
        kpoints
        incar keywords  npar, ispin, magmom, etc
            magmom will be list|dict such as ispin=2 & magmom=['N',2,...]|{'N':2,...(not checked)} 
    '''
    
    ### 1. Make atoms
    ## read poscar
    if len(poscar) == 1 and re.search('pos', poscar[0], re.I):
        atoms = io.read_poscar(poscar)
    ## generate surface
    elif len(poscar) >= 3:
        asize = int(poscar[2])
        size=(asize,asize,asize)
        print(f"fccsurf: {poscar[0]} {poscar[1]}, {size}")
        atoms = surflab.fccsurfaces(poscar[0], poscar[1], size, vac=15)
    else:
        atoms = None

    ### 2. Set params 
    if cat_kind == 'orr':
        ### INCAR params: 
        ###     magmom = dict or list: ispin=2, magmom=['N',2]
        incar_params = dict(npar=npar, kpoints=[1,1,1], ediff=0.0001, ediffg=-0.05, encut=400, ispin=2)
        sim_params   = dict(nproc=nproc)
        sim_params.update(incar_params)
    elif cat_kind == 'her':
        atoms = surflab.fccsurfaces('Pt', '111', (3,3,3), vac=15)
        sim_params  = dict(npar=npar, kpoints=[4,4,1], nproc=nproc, ediff=0.0001, ediffg=-0.04, encut=400)
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
            ### runORR: calc, sim_param, mode='opt', fix=None, active=None, vib=1, label='test'
            catalysis.runORR(calc, sim_params, mode=mode, vib=Lvib, fix='b1L', label=flabel)
        elif cat_kind == 'her':
            catalysis.runHER(calc, sim_params, mode=mode, vib=Lvib, label=flabel)
    elif job == 'incar':
        for k, v in calc.get_options():
            print(f"{k:>10}\t{v}")
    elif job == 'plot':
        ## To plot: get values of totE, zpe, TS and plot
        # 1. get DFT values
        totE, zpe, TS = catalysis.runORR(at, nproc=24, npar=4, mode='opt', kpoints=[4,4,1], vib=1, label='test')
        print(totE)
        print(zpe)
        print(TS)

        G_ORR      = catalysis.gibbs_ORR_4e_acid(TE=totE, pH=0)
        print(G_ORR)
        G_ORR_vib  = catalysis.gibbs_ORR_4e_acid(TE=totE, ZPE=zpe, TS=TS, pH=0)
        print(G_ORR_vib)
        catalysis.plot_ORR_4e_acid(G_ORR_vib, U=0.7, legend=['U=1.23V', 'U=0.70V', 'U=0.00V'])
    else:
        pass
    return 0

def main():
    parser = argparse.ArgumentParser(description="Running catalysis in slurm::\
                        \n\tselect job, subjob, some options for vib, overwrite\
                        \n\tsystem params is applied to this system")
    parser.add_argument('-j', '--job', default='run', choices=['run', 'plot', 'incar'], help='incar: show default params')
    parser.add_argument('-c', '--cat_kind', default='orr', choices=['orr', 'her', 'oer'], help='catalytic reactions')
    parser.add_argument('-l', '--flabel', default='test', help='label for dirname')
    parser.add_argument('-o', '--overwrite', action='store_true', help='if there exists dir, overwrite')
    parser.add_argument('-p', '--poscar', nargs='*', help="use any poscar or generate surface: ['Pt', '111', (3,3,3)]")
    group_vasp  = parser.add_argument_group(title='NC running options')
    group_vasp.add_argument('-m', '--mode', default='opt', choices=['opt', 'sp'], help='Opt mode')
    group_vasp.add_argument('-nv', '--novib', action='store_false', help='run orr without vibration')
    group_sys   = parser.add_argument_group(title='System-dependent inputs')
    group_sys.add_argument('-N', '--nnode', default=1, type=int, help='number of nodes: if needed')
    group_sys.add_argument('-np', '--nproc', type=int, default=24, help='number of process for mpirun')
    group_sys.add_argument('--npar', type=int, default=4, help='value in INCAR')
    parser.add_argument('-u', '--usage', action='store_true', help='explains how to run.')

    args = parser.parse_args()
    if args.usage:
        print(f"Usage::\
                \n\tThis is an example of a job submit in queue system such as slurm\
                \n\tRun sbatch with jobname, partition, nnode, nproc with variables\
                \n\t    sbatch -J test -p X3 -N 1 -n 20 --export=job='{args.cat_kind}' slurm_sbatch_nc.sh\
                \n\trun_catalysis.py is run inside job script\
                \n\t    run_catalysis.py -j orr -sj run -N {args.nnode} -np {args.nproc} --npar $npar\
                \n\tjob is running in work dir(jobname) & logfile is written in submit dir\
                \n\t    mpirun -np {args.nproc} VASP_EXC # in work dir\
                \n\tAfter job finished: jobid.jobname.log -> jobid.jobname.out\
            ")
        sys.exit(0)
    run_catalysis(args.job, args.cat_kind, args.flabel, args.overwrite, args.poscar, args.mode, args.novib, args.nnode, args.nproc, args.npar)

if __name__ == "__main__":
    main()
