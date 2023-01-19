import argparse
from nanocore import io
from nanocore import catalysis     

def run_catalysis(job, nnode, nproc, npar):
    '''
    Dir is overwritten
    If OUTCAR_job_label exists, cal for model is skipped
    Control params:
        mode    'opt', 'sp'(default) for optimization or single point calculation
        vib     activate vibration calculation for adsorbate for TS 
        label   filename of job to save OUTCAR, CONTCAR, XDATCAR, etc
        fix     
        active
    Vasp params for INCAR, KPOINTS & mpirun:
        nproc
        kpoints
        incar keywords  npar, ispin, magmom, etc
    '''

    if job == 'run':
        at = io.read_poscar('POSCAR')

        ### bind vasp kw to dict
        vasp_kw=dict(nproc=nproc, npar=npar, kpoints=[2,2,1] )

        ### Run VASP with analysis
        #catalysis.runORR(at, nproc=24, npar=4, mode='opt', kpoints=[4,4,1], vib=1, label='test')
        catalysis.runORR(at, vasp_kw, mode='opt', vib=1, label='test')
    elif job == 'incar':
        from nanocore.simulator.vasp import Vasp
        atoms       = io.read_poscar('POSCAR')
        for k, v in Vasp(atoms).get_options():
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
    return 0

def main():
    parser = argparse.ArgumentParser(description="Running catalysis in slurm")
    #parser.add_argument('-x', '--xpartition', help='if needed, specify nodename')
    parser.add_argument( '-j', '--job', default='run', choices=['run','plot','incar'], help='incar: show default params')
    parser.add_argument('-N', '--nnode', default=1, type=int, help='number of nodes: if needed')
    parser.add_argument('-np', '--nproc', type=int, default=24, help='number of process for mpirun')
    parser.add_argument( '--npar', type=int, default=4, help='value in INCAR')

    args = parser.parse_args()

    #run_catalysis(args.xpartition, args.nnode, args.nproc, args.npar)
    run_catalysis(args.job, args.nnode, args.nproc, args.npar)
    

if __name__ == "__main__":
    main()
