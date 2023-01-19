import argparse
from nanocore import io
from nanocore import catalysis     
from nanocore.simulator.vasp import Vasp

def print_params(**kw):
    for k, v in kw.items():
        print(f"{k:>10}\t{v}")
    return 0

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
    ### Make Vasp instance and pass to run_catalysis(Vasp_ins, cat='orr', 
    atoms   = io.read_poscar('POSCAR')
    calc    = Vasp(atoms)
    print_params(**dict(calc.get_options()))
    calc2    = Vasp(atoms)
    vas_param=dict(npar=4, kpoints=[2,2,1])
    calc2.set_options(**vas_param)
    print_params(**dict(calc2.get_options()))

    ### Run VASP with analysis
    if job == 'run':
        #catalysis.runORR(at, nproc=24, npar=4, mode='opt', kpoints=[4,4,1], vib=1, label='test')
        catalysis.runORR(calc2, vas_param, nproc=24,  mode='opt', vib=1, label='test')
    elif job == 'incar':
        for k, v in vas_ins.get_options():
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
    parser = argparse.ArgumentParser(description="Running catalysis in slurm")
    parser.add_argument( '-j', '--job', default='run', choices=['run','plot','incar','pass'], help='incar: show default params')
    parser.add_argument('-N', '--nnode', default=1, type=int, help='number of nodes: if needed')
    parser.add_argument('-np', '--nproc', type=int, default=24, help='number of process for mpirun')
    parser.add_argument( '--npar', type=int, default=4, help='value in INCAR')

    args = parser.parse_args()

    run_catalysis(args.job, args.nnode, args.nproc, args.npar)
    

if __name__ == "__main__":
    main()
