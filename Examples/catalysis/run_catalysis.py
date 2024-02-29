import argparse
import re, os, sys
from oorinano import catalysis
from oorinano.calculator.vasp import Vasp
from oorinano.calculator.vasp import readAtomicStructure as read_geo
from oorinano import surflab
from oorinano.utils.auxil import fname_ext, fname_root
import json

'''
input   Atomic structure should be given in a file
        metals are series for generating metal slabs
'''

def get_inputfile(inf, prefix):
    '''
    Get json input data file
    '''
    if inf:
        if os.path.exists(inf) and re.search('json', inf):
            infile = inf
        else:
            print("use json file for plot with -i ")
            sys.exit(10)
    else:
        cwd   = os.getcwd()
        files = os.listdir(cwd)
        jsonf = [ f for f in files if 'json' in f ]
        if len(jsonf) == 1 or len(jsonf) == 2:
            if len(jsonf) == 2 and prefix:
                fs = [ f for f in jsonf if prefix in f]
                infile = fs[0]
            else:
                infile = jsonf[0]
        else:
            print(f"There are abnormal number of json files in {cwd}: {jsonf}")
            sys.exit(11)
        return infile

def set_simulation_params(vasp_parallel, nproc):
    '''
    User defined parameters for INCAR, KPOINTS, and simulation such as nproc (number of process)
    User defined parameters will overwrite the default values the default values
    '''
    vasp_params={}
    ediff = 0.0001; ediffg=-0.05; encut=400
    ispin = 2
    kpoints=[4,4,1]
    ## npar vs ncore exclusive
    if re.match('c', vasp_parallel ):
        ncore=int(vasp_parallel[1:])
        vasp_params['ncore'] = ncore
    else:
        npar=int(vasp_parallel[1:])
        vasp_params['npar']  = npar
    
    sim_params = dict(nproc=nproc)

    ### orr or oer: vasp_params.update(dict(kpoints=[4,4,1], ediff=0.0001, ediffg=-0.05, encut=400, ispin=2))
    ###     magmom = dict or list: ispin=2, magmom=['N',2]
    #vasp_params['kband'] = "67 68 69 70"
    ### her: vasp_params.update(dict(kpoints=[4,4,1], ediff=0.0001, ediffg=-0.05, encut=400, ispin=2))
    ### add user defined parameters here
    vasp_params.update(dict(kpoints=kpoints, ediff=ediff, ediffg=ediffg, encut=encut, ispin=ispin))
    sim_params.update(vasp_params)
    return sim_params


def run_catalysis(cat_struct, surf_size, fix, job, cat_rxn, pH, flabel, nnode, nproc, Vasp_par ):
    '''
    cat_struct   file with any type for atomic structure | metal list
    job         run, model (check POSCAR), incar (check INCAR)
    cat_rxn     ORR(default=orr), HER, OER
    flabel      filename of job to save OUTCAR, CONTCAR, XDATCAR (default='test')
    
    Run params:
        flabel  tag for filename to save OUTCAR, CONTCAR, XDATCAR, etc
        fix     None (default)
                b1L for slab to fix bottom 1 layer of slab
        pivot   adsorbate anchored position: list (position) or (int) atom index
                default: 1. the highest z, 2. center of xy-plane of supercell
    VASP params for INCAR, KPOINTS & mpirun:
        nproc
        kpoints
        INCAR:  mode Lvib Vasp_par keywords  npar, ispin, magmom, etc
            magmom will be list|dict such as ispin=2 & magmom=['N',2,...]|{'N':2,...(not checked)}
            mod change value in Vasp.wrtie_INCAR, Vasp.run_calculator
            ispin define here (below)
    Other controllable parameters
        runORR()    calc, sim_param, mode, fix, pivot, vib, flabel, pH, cat_rnx=(orr|oer)
        runHER()    calc, sim_param, mode, fix, pivot, vib, flabel, pH, atoms_list=(many catalysts)
                    pivot  int-atom index, list (len=3): coordinate
                    oer     uses runORR

    '''
    ###### 1. Make AtomsSystem instances
    ### if cat_struct == str: it's structure file of any type
    if type(cat_struct) == str:
        inp = cat_struct
        struct_type = 'file'
        ### read POSCARs
        if re.search('pos', inp , re.I) or re.search('con', inp, re.I):
            atoms = read_geo(inp)
        # if different format add more to differentiate the file type
        else:
            print(f"input file error with {inp}")
        
    ### if cat_struct == list: generate metal slab
    else:
        struct_type = 'gen'
        liatoms=[]
        for metal in cat_struct:
            if len(surf_size) == 2:
                asize = int(surf_size[1])
                size=(asize,asize,asize)
            elif len(surf_size) == 4:
                size=(int(surf_size[1]), int(surf_size[2]), int(surf_size[3]))
            print(f"fccsurf: {metal} {surf_size[0]}, {size}")
            atoms = surflab.fccsurfaces(metal, surf_size[0], size, vac=15)     # bccsurfaces is also available 
            liatoms.append(atoms)
        atoms = liatoms[0]
        ### for slab structure, fix bottom layer
        if fix is None:
            fix = 'b1L'

    ###### 2. Set params: VASP params (INCAR, KPOINTS), Server params (nproc)
    sim_params = set_simulation_params(Vasp_par, nproc)

    ### Make Vasp instance and pass to runORR to make Vasp instance inside module 
    calc    = Vasp(atoms)
    calc.set_options(**sim_params)

    ### 3. Run catalysis (VASP) | Show INCAR | Plot
    ### to test set mode = 'sp', vib = False
    mode    = 'sp'; vib     = True                                                                 # user defined parameters
    if job == 'run':
        if cat_rxn == 'orr' or cat_rxn ==  'oer':
            catalysis.runORR(calc, sim_params, mode=mode, vib=vib, fix=fix, flabel=flabel, pH=pH, cat_rxn=cat_rxn)    #pivot = 24 (atom index)
        elif cat_rxn == 'her':
            if struct_type == 'file':
                catalysis.runHER(calc, sim_params, mode=mode, vib=vib, fix=fix, flabel=flabel, pH=pH)
            ### Make loop for several catalysts: make sub directories
            elif struct_type == 'gen':
                if 'liatoms' in locals():
                    catalysis.runHER(calc, sim_params, atoms_list=liatoms, mode=mode, fix=fix, flabel=cat_struct, pH=pH)
                else:
                    catalysis.runHER(calc, sim_params, mode=mode, fix=fix, flabel=flabel, pH=pH)
    elif job == 'model':
        calc.write_POSCAR()
    elif job == 'incar':
        for k, v in calc.get_options():
            print(f"{k:>10}\t{v}")
    return 0

def main():
    parser = argparse.ArgumentParser(description="Running catalysis::\
                        \n\tselect catalytic job, subjob [run, show incar, ...], some options for vib, overwrite\
                        \n\tsystem params partition, node, etc are applied to specific system")
    parser.add_argument('-j', '--job', default='run', choices=['run', 'model','incar', 'plot'], help='incar: show default params')
    parser.add_argument('-r', '--cat_rxn', default='orr', choices=['orr', 'her', 'oer'], help='catalytic reactions')
    parser.add_argument('-l', '--flabel', default='test', help='label for filename')
    struct = parser.add_mutually_exclusive_group()
    struct.add_argument('-i', '--inf', help="any type of atomic structure gets only single file")
    struct.add_argument('-g', '-m', '--gen', nargs='*', help='input metal series metal series')
    gplot = parser.add_argument_group(title='user-defined paramters: pH & U')
    gplot.add_argument('-ph', '--pH', default=0, type=int, help='get pH from command line')
    gplot.add_argument('-U', '--cellU', nargs='*', type=float, help='series of cell potential')
    gplot.add_argument('-c', '--colors', nargs='*', help = 'change color order')
    group_gen  = parser.add_argument_group(title='metal slab generation')
    group_gen.add_argument('-ss', '--surf_size', nargs='*', help=" ['111', 3, 3, 3]=[surface index & [sizex, sizey] size_z]")
    group_gen.add_argument('-fix', '--fix', help='fix bottom layer for slab')
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
                \n\t\t$python ../run_catalysis.py -j run -i POSCAR -np {args.nproc} --npar {args.npar} [--ncore $ncore]\
                \n\t\t$python ../run_catalysis.py -j run -r her -m Pt Au Ag Pd Ni Cu -ss 111 3 -np {args.nproc} --npar {args.npar} [--ncore $ncore]\
                \n\t3. Plot in the job directory\
                \n\t\t$python ../run_catalysis.py -r orr[oer] -ph 14\
            ")
        sys.exit(9)
    if args.ncore:
        nparallel='c'+str(args.ncore)
    else:
        nparallel='p'+str(args.npar)

    ### if read json file, plot here
    ### read json file for data (totE, zpe, TS) and replot with U and pH
    if args.job == 'plot':
        infile = get_inputfile(args.inf, args.cat_rxn)
        print(f"{infile}")
        with open(infile) as f:
            ene = json.load(f)
            if args.cat_rxn == 'orr' or args.cat_rxn == 'oer':
                if args.cat_rxn == 'orr':
                    gibbs_pH = catalysis.calc_gibbs_ORR_4e_pH(totE=ene['total energy'], zpe=ene['zpe'], TS=ene['TS'], pH=args.pH, Temp=298.15)
                elif args.cat_rxn == 'oer':
                    gibbs_pH = catalysis.calc_gibbs_OER_4e_pH(totE=ene['total energy'], zpe=ene['zpe'], TS=ene['TS'], pH=args.pH, Temp=298.15)
                catalysis.plot_OXR_4e_wU(gibbs_pH, flabel='plot', pH=args.pH, Temp=298.15, U=args.cellU, cat_rxn=args.cat_rxn, colors=args.colors)
            elif re.search('her', infile, re.I ):
                print("not coded for HER")
                sys.exit(3)
                #Gibbs_vib   = calc_gibbs_HER(totE_sys, totE_sysH, zpe, TS, pH=pH, Temp=T)
        return 0
   
    if args.inf:
        if os.path.exists(args.inf):
            atomic_struct = args.inf
        else:
            print(f"input file {args.inf} does not exist")
            sys.exit(2)
    elif args.gen:
        ### type(args.metals) == list
        atomic_struct = args.gen
    else:
        print(f"input atomic structures are given by -i POSCAR or -m metals for generation")
        sys.exit(1)
    ### input structure can be poscar (str) or generation metal (list)
    run_catalysis(atomic_struct, args.surf_size, args.fix, args.job, args.cat_rxn, args.pH, args.flabel, args.nnode, args.nproc, nparallel)
    return 0

if __name__ == "__main__":
    main()
