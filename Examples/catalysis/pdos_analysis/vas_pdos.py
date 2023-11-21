##!/home/joonho/anaconda3/bin/python
import oorinano as oori
import argparse
import sys, re
from oorinano.visualizer.mplottable import plot_dostable
from oorinano.calculator.vasp.file_format import *

def redefine_spin4file(fname):
    with open(fname, 'r') as f:
        ncol = len(f.readline().strip().split())
        if ncol == 5:
            Lspin = False
        elif ncol == 9:
            Lspin = True
        else:
            print(f"Wrong PDOS file format {fname}:: ncol {ncol} should be 5 or 9")
            sys.exit(21)
    return Lspin


def post_process(job, inf, latoms, Lspin, fname_pdos, orb, lplots, plot_option):
    '''
    job     keyword for job
        pd  pdos extraction
        an  orbital analysis such as band center, band maxdos, Fermi abundance
        pl  plot pdos and band center etc
    '''
    ###### File format: spin (un)polarized
    ### unploarized: iene s    p    d    t
    ### poloarized : iene s_up s_dn p_up p_dn d_up d_dn t_up t_dn
    

    ### Make ylegends to be plotted
    ### list of 0,1,2,...,s,s_up,s_dn,p,p_up,p_dn,d,d_up,d_dn,up,dn,all
    ### redefine Lspin using input file of PDOS

    
    


    if re.search('pd', job):    # pdos calculation
        if not latoms:
            print(f"to analyze PDOS, input atom indices using -a list")
            sys.exit(1)
        oori.vasp.calc_pdos(fname=inf, atom_list=latoms, Lspin=Lspin)

    if re.search('an', job):    # pdos analysis for dosmax, etc
        ### Anal: change orbital symbol to be analyzed into column indices
        ### For analysis: change into column index
        #if re.search('an', job):
        ### Redefine Lspin using PDOS file

        Lspin = redefine_spin4file(fname_pdos)
 
        if isinstance(orb, str) :   # if column, it is int
            ### For polarized cal
            if Lspin:
                if orb in pdos_polar_index.keys():
                    corb = pdos_polar_index['orb']
                else:
                    print(f"select one of orbital index for polarized cal: {pdos_polar_index}")
                    sys.exit(2)
            ### For unpolarized cal or sum of up, dn
            else:
                if orb.isdigit():
                    lorb = l2unpolar_orbital[orb]
                else:
                    lorb = orb
                corb = pdos_unpolar_index[lorb]
        else:
            corb = orb
        band_center, Ef_abund, Edosmax = oori.vasp.pdos_orbital_analysis( fname=fname_pdos, orb=corb)
        vert_line={'band center': band_center, 'dosmax': Edosmax[0]}
    else:
        vert_line={}
        lorb = None

    if re.search('pl', job):    # plot with control
        Lspin = redefine_spin4file(fname_pdos)
        ylegends=[]
        iys=[]
        if Lspin:
            if len(lplots) == 1:
                if lplots[0] == 'up':
                    ylegends = ['t_up', 's_up', 'p_up', 'd_up']
                    iys      = [  8,       2,       4,     6]
                    #Lsplit    = False
                elif lplots[0] == 'dn':
                    ylegends = ['t_dn', 's_dn', 'p_dn', 'd_dn']
                    iys      = [  9,       3,       5,     7]
                    #Lsplit   = False
                elif re.search('a', lplots[0]):
                    ylegends = [ 's_up', 'p_up', 'd_up', 's_dn', 'p_dn', 'd_dn']
                    iys      = [ 2,       4,     6,      3,       5,     7 ]
                    #Lsplit   = True
                else:
                    for porb in lplots:
                        ylegends.append(porb)
                        iys.append(pdos_unpolar_index[porb])
        else:
            if len(lplots) == 1 and re.search('a', lplots[0]): 
                ylegends = ['t', 's', 'p', 'd']
                iys      = [ 5,   2,   3,   4]
            else:
                for porb in lplots:
                    ylegends.append(porb)
                    iys.append(pdos_unpolar_index[porb])
    
        plot_dostable(fname_pdos, 1, iys,title='PDOS', xlabel=r'E-E$_F$ [eV]',ylabel='DOS',ylegend_in=ylegends,colors=None,lvertical=vert_line, plot_option=plot_option, orb=lorb, xlim=[-21, 15])
            
    return 0

def main():
    parser = argparse.ArgumentParser(description='Post process for VASP')
    parser.add_argument('-j', '--job', default='pdanpl', help='use pd for cal_pdos, an for orb analysis, pl for plot')
    parser.add_argument('-i', '--inf', default='DOSCAR', help='DOSCAR to be analyzed')
    parser.add_argument('-a', '--atom_list', type=int, nargs='*', help='atom list: ase-index + 1')
    parser.add_argument('-s', '--opt_spin', action='store_true', help='in case spin=2, sum or separate for up and down spin')
    glanal=parser.add_argument_group(title="ortibal analysis")
    glanal.add_argument('-pdos', '--pdoscar', default='SUM_ATOM.dat', help='fname for PDOS file')
    gorb=glanal.add_mutually_exclusive_group()
    gorb.add_argument('-l', '--anal_lorbital', help='l up or dn for analysis:[0|1|2|s|s_up|s_dn|p|p_up|p_dn|d|d_up|d_dn]')
    gorb.add_argument('-c', '--anal_column', type=int, help='column index of SUM_ATOM.dat' )
    gplot=parser.add_argument_group(title="select ortial to analysis and plot")
    gplot.add_argument('-p', '--plot_orbs', nargs='*', help='list of 0,1,2,...,s,s_up,s_dn,p,p_up,p_dn,d,d_up,d_dn,up,dn to be plotted')
    gplot.add_argument('-po', '--plot_option', default='split', help='plot up/dn split or opposite side')
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
            \n\t\t-l(c) angular momentum for analysis, only one orbital is analized\
            \n\t\t-s    spin option: in case of spin polarized cal, -s leaves up & down otherwise sum\
            \n\t\t-p    plot orbitals: list of 0,1,2,...,s,s_up,s_dn,p,p_up,p_dn,d,d_up,d_dn,up,dn to be plotted\
            \n\t\t-po   plot option for spin up & dn, 'sp'lit or drawn in 'op'posite side\
            \n\t    E.g.:\
            \n\t\t    python vas_pdos.py -a 43 -l p\
            \n\t\tpdos calculation\
            \n\t\t    python vas_pdos.py -j pd    -a 43\
            \n\t\torbital analysis\
            \n\t\t    python vas_pdos.py -j pdan  -a 43   -l p\
            \n\t\t    python vas_pdos.py -j an            -l 1\
            \n\t\tplot orbital\
            \n\t\t    python vas_pdos.py -j pdanpl -a 43  -l p -p a\
            \n\t\t    python vas_pdos.py -j anpl          -l p -p a\
            \n\t\t    python vas_pdos.py -j pl                 -p t s p\
            \n\t\t    python vas_pdos.py -j pdpl -a 43  -p a -s -po op[sp]\
            ")
        sys.exit(0)
    
    if args.anal_lorbital:
        orb = args.anal_lorbital    # '0' is letter
    elif args.anal_column:
        orb = args.anal_column      # 0 is integer
    else:
        orb = None
    
    
    post_process(args.job, args.inf, args.atom_list, args.opt_spin, args.pdoscar, orb, args.plot_orbs, args.plot_option)
    
    return 0


if __name__ == "__main__":
    main()