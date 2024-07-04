import argparse
import re, os, sys
from oorinano import catalysis
from oorinano.calculator.vasp import Vasp
from oorinano.calculator.vasp import readAtomicStructure
from oorinano import surflab
from oorinano.utils import np_Xn, host
import json

'''
input   Atomic structure should be given in a file
        metals are series for generating several metal slabs
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
    All the params defined here by USER
        incar_params    INCAR file
        vasp_params     KPOINTS, 
        sim_params      nproc, 

    vasp_parallel: c10|p12, VASP parallel key needs to be define
        c   NCORE   10
        p   NPAR    12
    nproc is passed to "mpirun vasp_executable -np num"

    User defined parameters for INCAR, KPOINTS, which will overwrite the default values
        ediff, ediffg, encut, ispin, ...
        ispin = 2   generates MAGMOM
    
    ** Special Jobs: U-correction, Solvent model, etc
    U-correction
        ldau    = True|False    only type 2 U-J is applied such as U-J 1.8  
        ldauuj  = value         ldauj = 1, ldauu = ldauj + value
    
    Solvent default: 78.3 for water
        lsol    = True
        EB_K    = value (optional)
    '''
    ### Overwrite INCAR default values
    ediff   = 0.0001; ediffg=-0.05; encut=400
    ispin   = 2
    incar_params=dict(ediff=ediff, ediffg=ediffg, encut=encut, ispin=ispin)
    
    ### Include Ucorr True|False
    ldau    = True 
    ldauuj  = 1.8    # this is for type 2 of LDAUU & LDAUJ
    if ldau:
        incar_params.update(dict(ldau=True, ldauuj=ldauuj))
    ### include solvent effect (water) or comment it out
    lsol = False     # python syntax
    if lsol:
        incar_params.update(dict(lsol=True))

    ### npar vs ncore exclusive
    if re.match('c', vasp_parallel ):
        ncore   = int(vasp_parallel[1:])
        incar_params['ncore'] = ncore
    else:
        npar=int(vasp_parallel[1:])
        incar_params['npar']  = npar
 
    vasp_params = dict(kpoints=[4,4,1])
    sim_params  = dict(nproc=nproc)

    sim_params.update(vasp_params)
    sim_params.update(incar_params)
    return sim_params


def run_catalysis(cat_structs, act_site, fix, job, cat_rxn, pH, flabel, nproc, Vasp_par):
    '''
    STRUCTURE input
        cat_struct  atoms list for 1 POSCAR or metal list for slab
        flabel      input or metal list
        act_site    atom index which starts from 1
        fix         option to fix layer in opt
    job         run, model (check POSCAR), incar (check INCAR)
    cat_rxn     ORR(default=orr), HER, OER
    flabel      filename of job to save OUTCAR, CONTCAR, XDATCAR (default='test')
    
    
    Run params:
        flabel  tag for filename to save OUTCAR, CONTCAR, XDATCAR, etc
        fix     None (default)
                b1L for slab to fix bottom 1 layer of slab
        act_site   adsorbate anchored position: list (position) or (int) which atom index start from 1
                default: 1. the highest z, 2. center of xy-plane of supercell
    VASP params for INCAR, KPOINTS & mpirun:
        nproc
        kpoints
        INCAR:  mode Lvib Vasp_par keywords  npar, ispin, magmom, etc
            magmom will be list|dict such as ispin=2 & magmom=['N',2,...]|{'N':2,...(not checked)}
            mod change value in Vasp.wrtie_INCAR, Vasp.run_calculator
            ispin define here (below)
        incar_arg   Ucorr, solvent
                    Ucorr: value
                    solvent: LSOL, EB_K=value (default: water=78.3) & SIGMA_K, NC_K 
    Other controllable parameters
        runORR()    calc, sim_param, mode, fix, pivot, vib, flabel, pH, cat_rnx=(orr|oer)
        runHER()    calc, sim_param, mode, fix, pivot, vib, flabel, pH, atoms_list=(many catalysts)
                    act_site  int-atom index which starts from 1, list (len=3): coordinate
                    oer     uses runORR

    '''
    ###### 1. Obtain prime atoms to get just info of model, incar
    prime_atoms = cat_structs[0]

    ###### 2. Get params: INCAR params, VASP params (KPOINTS), sim_params
    sim_params = set_simulation_params(Vasp_par, nproc)
    calc    = Vasp(prime_atoms)
    calc.set_options(**sim_params)
    
    ### to test set mode = 'sp', vib = False
    mode    = 'opt';    vib = True                                                                 # user defined parameters
    
    act_site = 40       # overwrite input: 19 for sac, 40 for apcc

    if job == 'run':
        if cat_rxn == 'orr' or cat_rxn ==  'oer':
            ### Make Vasp instance and pass to runORR to make Vasp instance inside module 
            catalysis.runORR(calc, sim_params, mode=mode, vib=vib, fix=fix, flabel=flabel, pH=pH, cat_rxn=cat_rxn, act_site=act_site)    #act_site = 24 (atom index)
        elif cat_rxn == 'her':
            catalysis.runHER(calc, sim_params, atoms_list=cat_structs, mode=mode, fix=fix, flabel=flabel, pH=pH)
    elif job == 'model':
        calc.write_POSCAR()
    elif job == 'incar':
        calc.run_calculator(mode=f"{mode}w")
    return 0

def main():
    parser = argparse.ArgumentParser(description="Running catalysis::\
                        \n\tselect catalytic job, subjob [run, show incar, ...], some options for vib, overwrite\
                        \n\tsystem params partition, node, etc are applied to specific system")
    parser.add_argument('-j', '--job', default='run', choices=['run', 'model','incar', 'plot'], help='incar: show default params')
    parser.add_argument('-r', '--cat_rxn', default='orr', choices=['orr', 'her', 'oer'], help='catalytic reactions')
    parser.add_argument('-l', '--flabel', default='test', help='label for filename')
    struct = parser.add_argument_group('structure-related args')
    struct.add_argument('-fix', '--fix', help="fix bottom layer for slab: 'b1L' for bottom 1 layer fixed")
    struct.add_argument('-as', '--active_site', type=int, help='atom index for active site which starts for 1')
    str_gen=struct.add_mutually_exclusive_group()
    str_gen.add_argument('-i', '--inf', help="any type of atomic structure gets only single file")
    str_gen.add_argument('-g', '-m', '--metals', nargs='*', help='input metal series metal series')
    struct.add_argument('-ss', '--surf_size', nargs='*', default = ['111', 3, 3, 3], help=" ['111', 3, 3, 3]=[surface index & [sizex, sizey] size_z]")
    gplot = parser.add_argument_group(title='user-defined paramters: pH & U')
    gplot.add_argument('-ph', '--pH', default=0, type=int, help='get pH from command line')
    gplot.add_argument('-U', '--cellU', nargs='*', type=float, help='series of cell potential')
    gplot.add_argument('-c', '--colors', nargs='*', help = 'change color order')
    group_sys   = parser.add_argument_group(title='System-dependent inputs')
    group_sys.add_argument('-q', '--qname', default='test', help='queue submit job name')
    group_sys.add_argument('-x', '--partition', default='X3', help='partition name as string')
    group_sys.add_argument('-n', '--nnode', default=1, type=int, help='number of nodes: if needed')
    group_sys.add_argument('-np', '--nproc', type=int, help='number of process for mpirun') # auto calculation from env_slurm
    parallel = group_sys.add_mutually_exclusive_group()
    parallel.add_argument('--npar', type=int, default=4, help='npar value in INCAR')
    parallel.add_argument('--ncore', type=int, help='ncore value in INCAR for KISTI')
    parser.add_argument('-u', '--usage', action='store_true', help='explains how to run.')

    args = parser.parse_args()
    #print(f"host {host} in {__name__}")
    if not args.nproc:
        if host == 'cluster':
            #print(f"np_Xn {np_Xn}")
            if 'np_Xn' in globals():
                nproc = args.nnode * np_Xn[args.partition]
            else:
                print("no args.nproc and no globals() for np_Xn")
        else:
            nproc = args.nnode * 40
    else:
        nproc = args.nproc
    
    if args.usage:
        print(f"Usage::\
                \n    These are examples for job submit in queue systems and direct run\
                \n    Check 'readme.txt' to set VASP envirionment\
                \n    Run:\
                \n\t1. Submit jobscript with jobname, (slurm: partition, nnode, nproc) with variables\
                \n\t    slurm::\
                \n\t\t$sbatch -J {args.qname} -p {args.partition} -N {args.nnode} -n {nproc}  slm_catalysis.sh\
                \n\t\t$sbatch -J {args.qname} -p {args.partition} -N {args.nnode} -n {nproc} --export=cat='{args.cat_rxn}',pos='gen' slm_catalysis.sh (default)\
                \n\t\t$sbatch -J {args.qname} -p {args.partition} -N {args.nnode} -n {nproc} --export=cat='{args.cat_rxn}',pos='POSCAR.xxx' slm_catalysis.sh\
                \n\t\t    - pos='POSCAR.xxx' to copy existing file to POSCAR & 'gen' for slab generation\
                \n\t    pbs  ::\
                \n\t\t$qsub -N {args.qname} pbs_vasp_kisti_skl.sh\
                \n\t\t$qsub -N {args.qname} -v cat='{args.cat_rxn}',pos='{args.inf}' pbs_vasp_kisti_skl.sh\
                \n\t    Output::\
                \n\t\t/test     job directory is generated\
                \n\t\trun_catalysis.py is run inside job script\
                \n\t\tjob is running in subdir(jobname) & logfile is written in workdir (submit dir)\
                \n\t\tjob finishes: jobname.log -> jobname.out\
                \n\t\tIn case rerun: missing OUTCAR is recalculated\
                \n\t2. Direct run inside job directory\
                \n\t    ORR with POSCAR\
                \n\t\t$python ../run_catalysis.py -j run -i POSCAR -np {args.nproc} --npar {args.npar} [--ncore $ncore]\
                \n\t\t$python ../run_catalysis.py -j run -r her -m Pt Au Ag Pd Ni Cu -ss 111 3 -np {args.nproc} --npar {args.npar} [--ncore $ncore]\
                \n\t3. Direct run to replot, print INCAR, check a generated POSCAR\
                \n\t\t$python ../run_catalysis.py -j run -r orr[oer] -ph 14\
                \n\t\t$python ../run_catalysis.py -j plot\
                \n\t\t$python ../run_catalysis.py -j incar -i POSCAR\
                \n\t\t$python ../run_catalysis.py -j model -m pt\
                \n\t\t$python ../run_catalysis.py -j model -m Pt  -ss 111 3\
            ")
        sys.exit(9)
    ### nparallel passes only ncore or npar
    if args.ncore:
        nparallel='c'+str(args.ncore)
    else:
        nparallel='p'+str(args.npar)

    ### 1. if job==plot, read json file
    ### read json file for data (totE, zpe, TS) and replot with U and pH
    if args.job == 'plot':
        infile = get_inputfile(args.inf, args.cat_rxn)
        print(f"{infile}")
        with open(infile) as f:
            ene = json.load(f)
            print(f"total energy: {ene['total energy']}")
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
    
    ###### 1. Make AtomsSystem instances
    lis_atoms=[]
    if args.inf:
        if os.path.exists(args.inf):
            atoms = readAtomicStructure(args.inf)
            lis_atoms.append(atoms)
        else:
            print(f"input file {args.inf} does not exist")
            sys.exit(2)
    ### metal list
    elif args.metals:
        args.flabel = args.metals       # metal name becomes file lable
        for metal in args.metals:
            if len(args.surf_size) == 2:
                asize = int(args.surf_size[1])
                size=(asize,asize,asize)
            elif len(args.surf_size) == 4:
                size=(int(args.surf_size[1]), int(args.surf_size[2]), int(args.surf_size[3]))
            print(f"fccsurf: {metal} {args.surf_size[0]}, {size}")
            atoms = surflab.fccsurfaces(metal, args.surf_size[0], size, vac=15)     # bccsurfaces is also available 
            lis_atoms.append(atoms)
        ### for slab structure, fix bottom layer
        if not args.fix:
            args.fix = 'b1L'
    
    ### input structure can be poscar (str) or generation metal (list)
    run_catalysis(lis_atoms, args.active_site, args.fix, args.job, args.cat_rxn, args.pH, args.flabel, nproc, nparallel)
    return 0

if __name__ == "__main__":
    main()
