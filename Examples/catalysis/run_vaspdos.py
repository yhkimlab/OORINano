import argparse, re, sys, os

from oorinano.calculator.vasp import Vasp, mod_vas
from oorinano.calculator.vasp import readAtomicStructure as read_geo
from oorinano.utils.auxil import fname_ext, fname_root
    
def cal_doscar(poscar, job, odir, nnode, nproc, Vasp_par):
    '''
    job         spchg   sp to generate CHGCAR
                dos     dir should have CHGCAR
                spdos   run two jobs (sp, dos) consecutively
    
    VASP params for INCAR, KPOINTS & mpirun:
        nproc
        kpoints
        poscar  Lvib
        INCAR:  mode Lvib Vasp_par keywords  npar, ispin, magmom, etc
            copy INCAR
            mod change value in Vasp.wrtie_INCAR, Vasp.run_calculator
    '''
    
    ### 1. POSCAR was copied in slm.sh
    atoms = read_geo(poscar)
    
    ### 2. INCAR read or make INCAR here
    old_incar = odir + '/INCAR'
    
    incar_params = mod_vas.read_incar(old_incar)
        
    sim_params   = dict(nproc=nproc)
    
    if re.match('c', Vasp_par ):
        ncore=int(Vasp_par [1:])
        incar_params['ncore'] = ncore
    else:
        npar=int(Vasp_par [1:])
        incar_params['npar']  = npar
        print(f"npar {incar_params['npar']}")
    
    calc    = Vasp(atoms)
     
    ### run | sp for CHGCAR | dos with CHGCAR | separately
    ### spdos: both can't be run at the same time because job script makes different directories for sp and dos
    if re.search('sp', job):
        if incar_params['ISPIN'] == 2:
            incar_params['MAGMOM'] = None   # MAGMOM would be different depending on POSCAR, if None, recalculat
        incar_params['lcharg'] = True     # Lwave = True
        sim_params.update(incar_params)
        calc.set_options(**sim_params)
        calc.run_calculator(mode='sp')
                
    elif re.search('dos', job):
        ### make directory in shell script
        incar_params.update(dict(kpoints=[8,8,1], nedos=4001, emin=-25, emax=15))            # increase kpoints for dos cal (post-process)
        sim_params.update(incar_params)
        calc.set_options(**sim_params)
        calc.run_calculator(mode='dos')
    else:
        print(f"Job error: select job in [sp, dos]")
        sys.exit(100)
    
    return 0

def main():
    parser = argparse.ArgumentParser(description="Running catalysis::\
                        \n\tselect catalytic job, subjob [run, show incar, ...], some options for vib, overwrite\
                        \n\tsystem params partition, node, etc are applied to specific system")
    parser.add_argument('-d', '-o', '--odir', help='old directory is used for slurm jobname')
    parser.add_argument('-p', '--poscar', default='CONTCAR_test_0_cat', help="select opt contcar in catalysis directory")
    parser.add_argument('-j', '--job', default='spdos', choices=['dos', 'sp', 'spdos'], help='sp for chgcar, dos for doscar')
    parser.add_argument('-u', '--usage', action='store_true', help='explains how to run.')

    group_sys   = parser.add_argument_group(title='System-dependent inputs')
    group_sys.add_argument('-x', '--partition', help='partition name')
    group_sys.add_argument('-n', '--nnode', default=1, type=int, help='number of nodes: if needed')
    group_sys.add_argument('-np', '--nproc', type=int, default=24, help='number of process for mpirun')

    parallel = group_sys.add_mutually_exclusive_group()
    parallel.add_argument('--npar', type=int, default=4, help='npar value in INCAR')
    parallel.add_argument('--ncore', type=int, help='ncore value in INCAR for KISTI')
    
    args = parser.parse_args()

    if args.usage:
        print(f"Usage::\
                \n    These are examples for job submit in queue systems and direct run\
                \n    Check 'readme.txt' to set VASP envirionment\
                \n    Run:\
                \n\t1. Submit jobscript with jobname(output dirname), (slurm: partition, nnode, nproc) with variables\
                \n\t    Postprocess: run dos\
                \n\t\t: Make CHGCAR and run dos\
                \n\t\t$sbatch -J {args.odir}{args.job} -p {args.partition} -N {args.nnode} -n {args.nproc} --export=job={args.job} slm_vaspdos.sh\
                \n\t\t$sbatch -J {args.odir}{args.job} -p {args.partition} -N {args.nnode} -n {args.nproc} --export=job={args.job},odir={args.odir},pos={args.poscar} slm_vaspdos.sh\
                \n\t\t: Only dos when aleardy CHGCAR exists\
                \n\t\t$sbatch -J {args.odir} -p {args.partition} -N {args.nnode} -n {args.nproc} --export=job='dos',odir={args.odir},pos={args.poscar} slm_vaspdos.sh\
                \n\t    pbs  ::\
                \n\t\t$qsub -N {args.odir}{args.job} pbs_vasp_kisti_skl.sh\
                \n\t    Output::\
                \n\t\tjob finishes: jobname.log -> jobname.out\
                \n\t2. Direct run inside job directory\
                \n\t\t$python ../run_vaspdos.py -j {args.job} -np 12 --npar 2\
                ")
        sys.exit(0)
    if args.ncore:
        nparallel='c'+str(args.ncore)
    else:
        nparallel='p'+str(args.npar)
            
    cal_doscar(args.poscar, args.job, args.odir, args.nnode, args.nproc, nparallel)

if __name__ == "__main__":
    main()
