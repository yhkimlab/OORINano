##!/home/joonho/anaconda3/bin/python
import oorinano as oori
import argparse
import sys, re
from oorinano.visualizer.mplottable import plot_dostable

def post_process(job, inf, latoms, option, fname_pdos, orb, iys, ylegends):
    '''
    job     keyword for job
        pd  pdos extraction
        an  orbital analysis such as band center, band maxdos, Fermi abundance
        pl  plot pdos and band center etc
    '''
    if re.search('pd', job):    # pdos calculation
        if not latoms:
            print(f"to analyze PDOS, input atom indices using -a list")
            sys.exit(1)
        oori.vasp.calc_pdos(fname=inf, atom_list=latoms, option=option)
    if re.search('an', job):    # pdos analysis for dosmax, etc
        band_center, Ef_abund, Edosmax = oori.vasp.pdos_orbital_analysis( fname=fname_pdos, orb=orb)
        vert_line={'band center': band_center, 'dosmax': Edosmax[0]}
    else:
        vert_line={}
    if re.search('pl', job):    # plot with control
        plot_dostable(fname_pdos, 1, iys,title='PDOS', xlabel=None,ylabel='DOS',ylegend_in=ylegends,colors=None,lvertical=vert_line, orb=orb, xlim=[-21, 15])
            
    return 0

def main():
    parser = argparse.ArgumentParser(description='Post process for VASP')
    parser.add_argument('-j', '--job', default='pdanpl', help='use pd for cal_pdos, an for orb analysis, pl for plot')
    parser.add_argument('-i', '--inf', default='DOSCAR', help='DOSCAR to be analyzed')
    parser.add_argument('-a', '--atom_list', type=int, nargs='*', help='atom list: ase-index + 1')
    parser.add_argument('-os', '--opt_spin', choices=['split', 'polar'], help='in case spin=2, sum, split or polar (- for down) for up and down spin')
    gorbial=parser.add_argument_group(title="ortibal analysis")
    gorbial.add_argument('-pdos', '--pdoscar', default='SUM_ATOM.dat', help='fname for PDOS file')
    gorb=gorbial.add_mutually_exclusive_group()
    gorb.add_argument('-l', '--lorbital', choices=['0', '1', '2', 's', 'p', 'd'], help='angular quantum number for analysis')
    gorb.add_argument('-c', '--column', help='column index of SUM_ATOM.dat' )
    gplot=parser.add_argument_group(title="select ortial to analysis and plot")
    gplot.add_argument('-po', '--plot_orb', default=['a'], nargs='*', choices=['s', 'p', 'd', 't', 'a'], help='use s,p,d,t for orbital to be plotted')
    #gplot.add_argument('-pl', '--dos_line', nargs='*', choices=['m', 'c', 'a', 'b'], help='options to plot, max, center, all')
    parser.add_argument('-u', '--usage', action='store_true', help='usage for main script')
    args = parser.parse_args()

    if args.usage:
        print(f"    Usage::\
            \n\tPrerequisit:    check atom index to plot pdos using POSCAR\
            \n\t                index start from 1\
            \n\t    Option:\
            \n\t\t-j    pdanpl  (default)\
            \n\t\t      pd      pdos calculation\
            \n\t\t      an      analysis such as dosmax, doscenter, fermi_abundance\
            \n\t\t      pl      plot option\
            \n\t\t-a    atom list\
            \n\t\t-l    angular momentum for analysis\
            \n\t\t-po   plot option: designate orbital to plot [default: 'a',all]\
            \n\t    E.g.:\
            \n\t\t    python vas_pdos.py -a 43 -l p\
            \n\t\tpdos calculation\
            \n\t\t    python vas_pdos.py -j pd    -a 43\
            \n\t\torbital analysis\
            \n\t\t    python vas_pdos.py -j pdan  -a 43   -l p\
            \n\t\t    python vas_pdos.py -j an            -l p\
            \n\t\tplot orbital\
            \n\t\t    python vas_pdos.py -j pdanpl -a 43  -l p -po a\
            \n\t\t    python vas_pdos.py -j anpl          -l p -po a\
            \n\t\t    python vas_pdos.py -j pl                 -po t s p\
            ")
        sys.exit(0)
    
    orbital = {'0': 's', '1': 'p', '2':'d'}
    ### orb = ['s', 'p', 'd', 'c#']
    if re.search('an', args.job):
        if args.lorbital:
            if args.lorbital.isdigit():
                orb = orbital[args.lorbital]
            else:
                orb = args.lorbital
        elif args.column:
            orb = 'c' + args.column
        else:
            print("input angular quantum number (-l) or column number (-c)")
            sys.exit(2)
    else:
        orb = None

    ### column and ylegend should match for spin sum file
    colm_lmatch={'s': 2, 'p': 3, 'd': 4, 't': 5}    # column index from 1 [x]
    ylegends=[]
    iys=[]
    if args.plot_orb[0] == 'a':
        ylegends = ['t', 's', 'p', 'd']
        iys      = [ 5,   2,   3,   4]
    else:
        for porb in args.plot_orb:
            ylegends.append(porb)
            iys.append(colm_lmatch[porb])

    post_process(args.job, args.inf, args.atom_list, args.opt_spin, args.pdoscar, orb, iys, ylegends)
    
    return 0


if __name__ == "__main__":
    main()