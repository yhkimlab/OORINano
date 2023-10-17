import oorinano as oori
import argparse

def post_process(job, inf, opt):
    if job == 'pdos':
        oori.vasp.pdos_orbital_analysis(fname=inf, orbitals=opt)
    elif job == 'split':
        at = oori.vasp.read_poscar(inf)       # in doscar_split.py
        n_at = len(at)
        oori.vasp.pdos_split_sum(sum_list=[i for i in range(n_at)])
    return 0

def main():
    parser = argparse.ArgumentParser(description='Post process for VASP')
    parser.add_argument('-j', '--job', choices=['pdos', 'split'], help='kind of post analysis')
    parser.add_argument('-i', '--inf', default='SUM_ATOM.dat', choices=['SUM_ATOM.dat','POSCAR'], help='input file')
    parser.add_argument('-o', '--option', default=1, type=int)
    args = parser.parse_args()

    pdos_analysis(args.job, args.inf, args.option)
    return 0


if __name__ == "__main__":
    main()

'''
doscar_split.py
at = oori.vasp.read_poscar('POSCAR')
n_at = len(at)
oori.vasp.pdos_split_sum(sum_list=[i for i in range(n_at)])
'''
