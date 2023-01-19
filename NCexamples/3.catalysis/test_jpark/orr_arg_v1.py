import argparse
from nanocore import io
from nanocore import catalysis     

def run_catalysis(xpart, nnode, nproc, npar):

    at = io.read_poscar('POSCAR')
    if not nproc:
        nproc=24
    if not npar:
        npar=4

    ### How to put in vasp keyword 
    kv_dict=dict(nproc=nproc, npar=npar, mode='opt', kpoints=[2,2,1], vib=1, label='test')

    ### Run VASP with analysis
    #catalysis.runORR(at, nproc=24, npar=4, mode='opt', kpoints=[4,4,1], vib=1, label='test')
    #catalysis.runORR(at, nproc=nproc, npar=npar, mode='opt', kpoints=[2,2,1], vib=1, label='test' )
    catalysis.runORR(at, kv_dict**)

    ### belows for explicit plot
    ### if each image is calculated in the separate directory, the calculation can be skipped if already calculation was done
    #TE, ZPE, TS = catalysis.runORR(at, nproc=24, npar=4, mode='opt', kpoints=[4,4,1], vib=1, label='test')
    #print(TE)
    #print(ZPE)
    #print(TS)

    #G_ORR      = catalysis.gibbs_ORR_4e_acid(TE=TE, pH=0)
    #print(G_ORR)
    #G_ORR_vib  = catalysis.gibbs_ORR_4e_acid(TE=TE, ZPE=ZPE, TS=TS, pH=0)
    #print(G_ORR_vib)
    #catalysis.plot_ORR_4e_acid(G_ORR_vib, U=0.7, legend=['U=1.23V', 'U=0.70V', 'U=0.00V'])
    return 0

def main():
    parser = argparse.ArgumentParser(description="Running catalysis in slurm")
    parser.add_argument('-x', '--xpartition', help='if needed, specify nodename')
    parser.add_argument('-N', '--nnode', default=1, type=int, help='number of nodes: if needed')
    parser.add_argument('-np', '--nproc', type=int, help='number of process: if needed')
    parser.add_argument( '--npar', help='value in INCAR')

    args = parser.parse_args()

    run_catalysis(args.xpartition, args.nnode, args.nproc, args.npar)
    

if __name__ == "__main__":
    main()
