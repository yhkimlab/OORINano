##!/home/joonho/anaconda3/bin/python
import oorinano as oori
import argparse
import sys

def post_process(job, inf, latoms, option, fname_pdos, orb, opt_plot):
    
    if job == 'pdos':
        oori.vasp.calc_pdos(fname=inf, atom_list=latoms, option=option)
        oori.vasp.pdos_orbital_analysis( fname=fname_pdos, orb=orb,plot_option=opt_plot)
    elif job == 'anal':
        oori.vasp.pdos_orbital_analysis( fname=fname_pdos, orb=orb,plot_option=opt_plot)
        ### plot
        
    return 0

def main():
    parser = argparse.ArgumentParser(description='Post process for VASP')
    parser.add_argument('-j', '--job', default='pdos', choices=['pdos', 'anal'], help='kind of post analysis')
    parser.add_argument('-i', '--inf', default='DOSCAR', help='DOSCAR to be analyzed')
    parser.add_argument('-a', '--atom_list', default=[1], type=int, nargs='*', help='atom list: ase-index + 1')
    parser.add_argument('-os', '--opt_spin', choices=['split', 'polar'], help='in case spin=2, sum, split or polar (- for down) for up and down spin')
    gorbial=parser.add_argument_group(title="ortibal analysis")
    gorbial.add_argument('-p', '--pdoscar', default='SUM_ATOM.dat', help='fname for PDOS file')
    gorb=gorbial.add_mutually_exclusive_group()
    gorb.add_argument('-l', '--lorbital', choices=['0', '1', '2', 's', 'p', 'd'], help='angular quantum number for analysis')
    gorb.add_argument('-c', '--column', help='column index of SUM_ATOM.dat' )
    gorbial.add_argument('-op', '--opt_plot', choices=['m', 'c', 'a', 'b'], help='options to plot, max, center, all')
    parser.add_argument('-u', '--usage', action='store_true', help='usage for main script')
    args = parser.parse_args()

    if args.usage:
        print(f"    Usage::\
            \n\tPrerequisit:    check atom list to plot pdos\
            \n\t                index start from 1\
            \n\tAbtain pdos:\
            \n\t    E.g.:\
            \n\t\tpython pdoscar.py -a 43\
            \n\t\tpython pdoscar.py -j plods -i DOSCAR -a atomlist -o split\
            \n\t    Options for pdos\
            \n\t\t-j  pdos    s, p, d sum for atom list\
            \n\t\t-i  'DOSCAR' (default)\
            \n\t\t-a  atom list   index starts from 1\
            \n\t\t-os spin option for spin = 2\
            \n\t\t    None    up and down spin summed: format(ene, s, p, d, s+p+d)\
            \n\t\t    split   separate up/down spin: format(ene, s_up, s_down, ... 9 colm)\
            \n\t\t    polar   down spin has - value\
            \n\tAnalysis of PDOS\
            \n\t    E.g.:\
            \n\t\tpdoscar.py -j anal -l 1\
            \n\t    Options for analysis\
            \n\t\t-j  anal    to get l-center, DOSmax, Fermi abundance\
            \n\t\t-p  PDOS file name 'SUM_ATOM.dat' (default)\
            \n\t\t-l  angular quantum number for analysis\
            \n\t\t-oa options for analysis\
            \n\t\t    max     energy for max DOS of l-orbital\
            \n\t\t    center  average E of l-orbital\
            \n\t\t    fa      E_fermi abundance\
            ")
        sys.exit(0)
    
    orbital = {'0': 's', '1': 'p', '2':'d'}
    ### orb = ['s', 'p', 'd', 'c#']
    if args.job == 'anal':
        if args.lorbital:
            if args.lorbital.isdigit():
                orb = orbital[args.lorbital]
            else:
                orb = args.lorbital
        elif args.column:
            orb = 'c' + args.column
        else:
            print("input angular quantum number (-l) or column number (-c)")
            sys.exit(1)


    post_process(args.job, args.inf, args.atom_list, args.opt_spin, args.pdoscar, orb, args.opt_plot)
    return 0


if __name__ == "__main__":
    main()

'''
doscar_split.py
at = oori.vasp.read_poscar('POSCAR')
n_at = len(at)
oori.vasp.pdos_split_sum(sum_list=[i for i in range(n_at)])
'''
