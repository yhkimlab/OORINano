import argparse
import os
from nanocore import io
from nanocore import catalysis     
from nanocore.simulator.vasp import Vasp
from nanocore import surflab

def print_params(**kw):
    for k, v in kw.items():
        print(f"{k:>10}\t{v}")
    return 0

def run_catalysis(job, subjob,  nnode, nproc, npar):
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
    ### Make Vasp instance and pass to runORR 
    if job == 'orr' and os.path.isfile('POSCAR'):
        atoms       = io.read_poscar('POSCAR')
        ### default: nproc=24, npar=4, mode='opt', kpoints=[4,4,1], vib=1, label='test'
        sim_params  = dict(npar=npar, kpoints=[2,2,1], nproc=nproc, ediff=0.0001, ediffg=-0.05, encut=400)
    elif job == 'her':
        atoms       = surflab.fccsurfaces('Pt', '111', (3,3,3), vac=15)
        ### default: npar=8, nproc=48, kpoints=[4,4,1], ediffg=-0.04
        sim_params  = dict(npar=npar, kpoints=[2,2,1], nproc=nproc, ediff=0.0001, ediffg=-0.04, encut=400)
    calc    = Vasp(atoms)
    calc.set_options(**sim_params)
    
    ### Run VASP with analysis
    if subjob == 'run':
        if job == 'orr':
            #catalysis.runORR(at, nproc=24, npar=4, mode='opt', kpoints=[4,4,1], vib=1, label='test')
            catalysis.runORR(calc, sim_params, mode='opt', vib=1, label='test')
        elif job == 'her':
            catalysis.runHER(calc, sim_params, mode='opt', vib=1, label='test')
    elif subjob == 'incar':
        for k, v in vas_ins.get_options():
            print(f"{k:>10}\t{v}")
    elif subjob == 'plot':
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
    parser.add_argument( '-j', '--job', default='orr', choices=['orr','her','oer'], help='kind of catalysis')
    parser.add_argument( '-sj', '--subjob', default='run', choices=['run','plot','incar'], help='incar: show default params')
    parser.add_argument('-N', '--nnode', default=1, type=int, help='number of nodes: if needed')
    parser.add_argument('-np', '--nproc', type=int, default=24, help='number of process for mpirun')
    parser.add_argument( '--npar', type=int, default=4, help='value in INCAR')

    args = parser.parse_args()

    run_catalysis(args.job, args.subjob, args.nnode, args.nproc, args.npar)
    

if __name__ == "__main__":
    main()
