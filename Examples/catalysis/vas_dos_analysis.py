#!/home/joonho/anaconda3/bin/python
import argparse
import sys, re
import oorinano as oori
from oorinano.visualizer.mplottable import plot_dostable

### output DOS data file format
ff_unpolar  = ['iene', 's', 'p', 'd', 'sum']    # sum = s + p + d
ff_polarsum = ff_unpolar
ff_polar    = ['iene', 's_up', 's_dn','p_up', 'p_dn', 'd_up', 'd_dn', 'sum_up', 'sum_dn']

def dos_analysis(job, doscar, latoms, Lspin, fname_pdos, an_orb, elimit, plot_spins, plot_style, ylegends):
    '''
    job     keyword for job
        pd  pdos extraction
        an  orbital analysis such as band center, band maxdos, Fermi abundance
        pl  plot pdos and band center etc
    '''
    
    ### 1. Lspin for un/polarized pdos calculation
    ### spin : autocalculated in module
    if re.search('pd', job):    # pdos calculation
        if not latoms:
            print(f"to analyze PDOS, input atom indices using -a list")
            sys.exit(1)
        ### returns Atom...dat, copy it to SUM_ATOM.dat, and make TDOS.dat
        oori.vasp.calc_pdos(fname=doscar, atom_list=latoms, Lspin=Lspin) 

    ### 2. PDOS analysis for dosmax, band center, Fermi abundance level
    ### If analysis for up and down, change plot_spins
    if re.search('an', job):    
        if an_orb in ff_unpolar:
            icol = ff_unpolar.index(an_orb)
        elif an_orb in ff_polar:
            icol = ff_polar.index(an_orb)
            if re.search('up', an_orb):         # plot up or dn
                plot_spins = 'up'
            elif re.search('dn', an_orb):
                plot_spins = 'dn'
        else:
            print(f"input orbital labe error {an_orb}")
        band_center, Ef_abund, Edosmax = oori.vasp.pdos_orbital_analysis(fname=fname_pdos, icol=icol, elimit=elimit)
        vert_line={'center': band_center, 'max': Edosmax[0], 'D_F': Ef_abund}    # key used as legend
    else:
        vert_line={}
        icol = None
    
    ### 3. Plot PDOS
    if re.search('pl', job):
        if Lspin:
            iys_up  = [ 7,  1,  3,  5]
            iys_dn  = [ 8,  2,  4,  6]
            if re.search('[ab]', plot_spins):   # 'all' or 'both
                iys = iys_up[1:] + iys_dn[1:]
            elif re.search('up', plot_spins):
                iys = iys_up
                #Lsplit    = False
            elif re.search('dn', plot_spins):
                iys = iys_dn
        else:
            iys     = [ 4,  1,  2,  3]
        print(f"iys {iys}")
        ### Make ylegends to be plotted
        if not ylegends:
            ylegends=[]
            for i in iys:
                if Lspin:
                    ylegends.append(ff_polar[i])
                else:
                    ylegends.append(ff_unpolar[i])
        print(f"plot {ylegends}")

        plot_dostable(fname_pdos, 0, iys,title='PDOS', xlabel=r'E-E$_F$ [eV]',ylabel='DOS',ylegends=ylegends,colors=None,lvertical=vert_line, plot_style=plot_style, orb=an_orb, xlim=[-25, 20])
            
    return 0

def main():
    parser = argparse.ArgumentParser(description='Post process for VASP')
    parser.add_argument('-j', '--job', default='pdanpl', help='pdos, analysis, and plot: calculate PDOS, orb analysis for dosmax, bandcenter, EFabundance, pl for plot')
    parser.add_argument('-i', '--infile', default='DOSCAR', help='DOSCAR to be analyzed')
    parser.add_argument('-a', '--atom_list', type=int, nargs='*', help='atom list: start from 1, cf. ASE from 0')
    parser.add_argument('-s', '--opt_spin', action='store_true', help='For spin [un]polarized data and plot')
    ganal=parser.add_argument_group(title="ortibal analysis")
    ganal.add_argument('-idos', '--pdosfile', default='SUM_ATOM.dat', help='fname for PDOS file')
    ganal.add_argument('-l', '--lorbital', default='d', help='combination of l and up/dn with format of l_up or l_dn ')
    ganal.add_argument('-el', '--elimit', default='inf', choices=['inf','fermi'], help='upper limit for integral')
    gplot=parser.add_argument_group(title="select ortial to analysis and plot")
    gplot.add_argument('-p', '--plot_spins', default='all', help='abupdn: all[both], up, dn')
    gplot.add_argument('-ps', '--plot_style', default=None, choices=['split', 'polar'], help='plot up/dn split or opposite side')
    gplot.add_argument('-yl', '--ylegends', help='can modify legends')
    parser.add_argument('-u', '--usage', action='store_true', help='usage for main script')
    args = parser.parse_args()

    if args.usage:
        print(f"    Usage::\
            \n\tPrerequisit:    check atom index to plot pdos using POSCAR\
            \n\t                index start from 1\
            \n\t    Option:\
            \n\t\t-j    pdanpl  (default)\
            \n\t\t      pd      pdos calculation\
            \n\t\t      an      analysis such as dosmax, doscenter, fermi_abundance level\
            \n\t\t      pl      plot option\
            \n\t\t-a    atom list\
            \n\t\t-s    spin option: makes PDOS with different filename\
            \n\t\t-l    angular momentum for analysis, only one orbital is analized\
            \n\t\t-p    plot orbitals: up|dn|both(all)\
            \n\t\t-ps   plot style for spin up & dn, 'sp'lit or drawn in 'polar' opposite side\
            \n\t    E.g.:\
            \n\t\t Non spin-polarized or summed for up- and dn-spin for polarized calculation without -s\
            \n\t\t    python ../vas_pdos.py -a 24\
            \n\t\t    python ../vas_pdos.py -a 43 -l p  # default = d\
            \n\t\t Spin-polarized: only one spin\
            \n\t\t    python ../vas_pdos.py -a 24 -s -l d_up\
            \n\t\t Spin-polarized: both spins without orbital analysis\
            \n\t\t    python ../vas_pdos.py -j pdpl -a 24 -s -p a -ps polar\
            \n\t\t    python ../vas_pdos.py -j pdpl -a 24 -s -p a -ps split\
            ")
        sys.exit(0)
    
    dos_analysis(args.job, args.infile, args.atom_list, args.opt_spin, args.pdosfile, args.lorbital, args.elimit, args.plot_spins, args.plot_style, args.ylegends)
    
    return 0


if __name__ == "__main__":
    main()
