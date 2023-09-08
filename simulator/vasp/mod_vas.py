### python-vasp module
""" 
    repository for many vasp script 
    get_vasp_repository()
    get_atoms_4pos() 
    make_mag_4pos()
    def make_incar(dic, iofile):
    def make_kpoints(kp, MH|gamma, **args=[band])
"""

import os
import json
import re
import sys
import argparse
#from common import whereami
from .. import atoms as ncatoms

### VASP ini & outfiles
vasf_default=['CHG','CHGCAR','CONTCAR','DOSCAR','EIGENVAL','IBZKPT','OSZICAR','OUTCAR','PCDAT','REPORT','vasprun.xml','WAVECAR',  'XDATCAR']
vasf_ini=['POSCAR','KPOINTS','INCAR','POTCAR']

### VASP jobs group with file to be modified
jg_poscar=['ini', 'zpe']      # ini uses, zpe modifies, others use CONTCAR
jg_kpoints=['dos','band','kp'] # kp for change kp file
jg_incar=['sp','opt','copt','vdw','chg','chgw','dos','pchg','band','mag','kisti']
jg_potcar=['lda','gga']
jg_link=['dos','band','pchg']

eps_H2O = 78.3
ini_dvasp = '/tmp'

TMUcorr = {'Fe': 1.8, 'Co': 2.1, 'Ni': 2.0, 'Pt':1.2}
natom = 999
incar_kw={  'ISTART': '0', 'ICHARG': '2', 'ISPIN': '1', 'MAGMOM': f'{natom} * 0',
            'GGA'   : 'PE', 'PREC': 'Accurate', 'ALGO': 'Fast', 'NPAR': 4, 'LPLANE': '.TRUE.',
            'ENCUT' : '400', 'NELMIN': '4', 'NELM' : '500', 'EDIFF': '1E-4', 'ISYM': '0', 'ADDGRID':'.TRUE.', 'LREAL':'Auto'
            }

def get_hostname():
    
    hostname = os.popen('hostname').read().rstrip()
    print(hostname)
    if hostname == 'login':
        hname = 'mlet'
    elif re.match('login0', hostname):
        hname = 'kisti'
    elif re.match('tgm-master',hostname):
        hname = 'pt'
    else:
        hname = hostname
    return hname    

def get_vasp_repository():
    """ 
        my vasp repository for POTCAR & KPOTINS 
    """
    global ini_dvasp

    ### subprocess is not working
    #hostname = subprocess.popen('hostname', stdout=subprocess.PIPE, shell=True)
    #proc = subprocess.Popen(["cat", "/etc/services"], stdout=subprocess.PIPE, shell=True)
    #(out, err) = proc.communicate()
    #print "program output:", out
    hostname = get_hostname()
    user=os.getenv('USER')
    #if hostname == 'chi' or hostname == 'pt' or hostname == 'iron' or hostname == 'mlet':
    if hostname == 'kisti':
        ini_dvasp = f"/home01/{user}/sandbox/pyvasp/ini"
    else:
        ini_dvasp = '/home/joonho/sandbox/pyvasp/ini'

    print(f"{__name__}:{whereami()}:: vasp repository is {ini_dvasp} in system {hostname}")
    if not os.access(ini_dvasp, os.F_OK):
        print("Error:: the directory cannot be found\n stop")
        exit(1)
    return ini_dvasp

def get_atoms_4pos(pos='POSCAR'):
    with open(pos, 'r') as f:
        lines = f.readlines()
        for index, line in enumerate(lines):
            if index <= 4: continue
            if line.strip().replace(' ','').isalpha():
                atoms=line.strip().split()
                continue
            if 'atoms' in locals():
                natoms=line.strip().split()
                break
    if 'atoms' in locals() and 'natoms' in locals():
        return natoms, atoms
    else:
        return 'err' 
                

def make_magmom_4pos(pos='POSCAR', magin=None):
    '''
    Read POSCAR [and input magmom]
    Return string of MAGMOM
    
    pos     POSCAR
    magin   MAGMOM input list of [ 'atom symbol', magmom, ... ]

    '''
    magmom = {}
    if magin:
        li = iter(magin)
        magmom = dict(zip(li, int(li)))
    natoms, atoms = get_atoms_4pos(pos)
    magstr="MAGMOM = "
    Lmag = False
    for index, atom in enumerate(atoms):
        if magin and atom in magin.keys():
            magstr += f"{natoms[index]}*{magin[atom]*1.5} "
            Lmag = True
        elif atom in ncatoms.atom_prop.keys():
            magstr += f"{natoms[index]}*{ncatoms.atom_prop[atom][1]*1.5} "
            Lmag = True
        else:
            print("ERROR: no magmom in input and repo")
            sys.exit(10)
            #magstr += f"{natoms[index]}*0 "
    if Lmag: return magstr + "100*0"
    else:    return "# " + magstr


def make_kpoints(kp, method):
    """ 
        Make KPOINTS file 
        only Gamma w. 1 1 1 and MH are adapted
    """
    fname = 'KPOINTS'

    if not kp:
        if method == 'gamma':
            kfile = ini_dvasp + '/kp.gamma'
            s = 'cp %s KPOINTS' % kfile
            print(f"{__name__}:{whereami()}::KPOINTS was copied from {kfile}")
            os.system(s)
            return 0
        else:
            print('more info for KPOINTS')

    f = open(fname, 'w')
    f.write("Automatic Mesh\n")                 # 1st line, description
    f.write("0\n")                              # 2nd line, number of K, 0 for automatic
    
    if re.match("g", method, re.IGNORECASE):    # 3rd line, Auto - gamma-centered MH Pack
        f.write("Gamma\n")                      #   gamma 
    else:
        f.write("Monkhosrt\n")                  #   MH
    if len(kp) == 3:                            # kp input as "1 1 1"
        l = "  ".join(kp) + "\n"                # 4th liine
        f.write(l)                              #
    else:
        kp_l = kp + "\n"
        f.write(kp_l)

    add = "000"                                 # 5th for shift of kp's
    l = "  ".join(add) + "\n"
    f.write(l)
    f.close()        
    return 0            

vasp_gga={'pbe': 'pe', 'rpbe': 'rp',  'revpbe': 're'}

def get_Ucorr(atoms):
    '''
    LDAU = .TRUE.
    LDAUTYPE = 2                     ! type1=LSDA+U (U & J), type2= LSDA+U (U-J)
    LDAUL  =  -1  -1   2  -1  -1  -1 ! No_corr(-1), p(1), d(2), f(3)
    LDAUU  =  0.0 0.0 2.8 0.0 0.0 0.0
    LDAUJ  =  0.0 0.0 1.0 0.0 0.0 0.0
    LDAUPRINT = 2        ! output occupation matrix.
    '''
    ldaul = 'LDAUL  = '
    ldauu = 'LDAUU  = '
    ldauj = 'LDAUJ  = '
    for atom in atoms:
        if atom in TMUcorr:
            ldaul += f'{2:^4}'
            ldauu += f'{TMUcorr[atom]+1.0:4.1f}'
            ldauj += f'{1.0:4.1f}'
        else:
            ldaul += f'{-1:^4d}'
            ldauu += f'{0.0:4.1f}'
            ldauj += f'{0.0:4.1f}'
    return ldaul, ldauu, ldauj
        
        

def make_incar(dic, iofile):
    """ Make INCAR file from dictionary of command line arguments
        automatically write "incar.key"-json string
        when modify "incar.key json string" 
        filename can be change with -f, --iofile """

    ### args list:  1. nonmag, mag(fm), afm
    ###             2. ini, cont
    ###             3. precision
    ###             4. gga & parallel
    ###             5. dispersion
    ###             6.

    ### if iofile exists, read, otherwise, write
    if os.path.isfile(iofile):
        dic = json.load(open(iofile))
        print(f"read {iofile} to make INCAR")
    else:
        with open(iofile, 'w') as ofile:
            ofile.write(json.dumps(dic, indent=4))
            print(f"write {iofile}; check it and run 'vas_make_incar.py' again")
            #sys.exit(1)

    print(dic)
    ################# PRE DEFINE ###############################
    fname = 'INCAR'
    f = open(fname, 'w')
    L_hybrid = 0
    L_vdw = 0
    L_paw = 0
    ########### DEFINE GGA name & hybrid
    if len(dic['dft'])==2:
        gga = dic['dft']
    elif len(dic['dft'])>2:
        gga = dic['dft'][:2]
        odd = dic['dft'][2:]
        if re.search('vdw', odd):
            L_vdw = 1
        if re.search('0', odd):
            L_hybrid=1
            dic['ini']='wav'
    print("hybrid %d vdW-DF %d" % (L_hybrid, L_vdw))
    ########### ENcut vs Precision
    L_encut=0
    if 'cutoff' in dic.keys():
        L_encut=1
    ########## WRITE option && MD option && SYSTEM options
    laechg = '.FALSE.'
    lwave = '.FALSE.'
    lcharg = '.FALSE.'
    ilog = dic['log']
    L_md = 0
    if dic['dynamics']:
        ilog = 0
        L_md = 1
        isym = 0
        kp = 'gamma'
        nkred = 0
        dic['precision']='normal'
        TK = 353
        tebeg = TK
        teend = TK
    ### system character
    elif dic['system'] == 'mol':
        isym = 0
        kp = 'gamma'
        nkred = 0
    else:
        nkred = 2
        #lreal = '.FALSE.'

    ########################## WRITE ############################################
    ### 0: system name
    comm = 'SYSTEM = ' + json.dumps(dic) + '\n\n'
    f.write(comm)
    ### 1: magnetism
    if dic['mag'] == 'nm':
        comm = '#MAGMOM:: natom*magnetism ... \n\n'
        f.write(comm)
    else:
        
        magmom = dic['magmom']
        comm = 'MAGMOM = '
        if dic['mag'] == 'fm':
            comm += str(natom) +'*'+magmom+' 999*0\n\n'
        else:
            for i in range(0, natom):
                if i % 2 == 0:
                    comm += str(magmom)+' '
                else:
                    comm += '-'+str(magmom)+' '
            comm += '999*0\n\n'                        
        f.write(comm)
    ### 2: start
    f.write('# continuation\n')
    if dic['ini'] == 'start' or L_md:
        istart = 0
        icharg = 2
    elif dic['ini'] == 'chg':
        istart = 1
        icharg = 1
    else:       # cont goes to wav
        istart = 1
        icharg = 0
    comm = 'ISTART = %d\nICHARG = %d\n' % (istart, icharg)    
            
    if dic['mag'] == 'nm':
        comm += 'ISPIN = 1\n\n'
    else:
        comm += 'ISPIN = 2\n\n'
    f.write(comm)

    ### 3: precision 
    f.write('# precision \n')
    if L_encut:
        com1 = 'ENCUT = %d\n' % dic['cutoff']
    if L_md:
        com1 += 'PREC = %s\n' % dic['precision']
    com1 += 'ISMEAR = 0 ; SIGMA = 0.05\n'
    com1 += 'NELMIN = 4 #; NELM = 500       # increase NELMIN to 4 ~ 8 in case MD|Ionic relax\n'
    if dic['relax'] == 'sp' or L_md:
        com1 += 'EDIFF = 1E-5  #; EDIFFG = -0.025 for MD or SP calc\n\n'
    else:
        com1 += 'EDIFF = 1E-5 ; EDIFFG = -0.025\n\n'
    com1 += 'ADDGRID = .TRUE.\n\n'
    f.write(com1)
    f.write('# mixing\n')
    com1 = 'MAXMIX = 40 for MD or ionic relaxation: reuse\n'
    com1 += '#IMIX = 4; #AMIX = 0.2; #BMIX = 0.0001; #AMIX_MAG = 0.8; #BMIX_MAG = 0.0001\n\n'
    if nkred:
        com1 += 'NKRED = %d\n' % nkred
    f.write(com1)
    f.write('# parallel performance and gga\n')
    #com1 = 'ALGO = Fast\n'
    #com1 += '#IALGO=48\n'
    #com1 += 'NSIM = 4; NPAR = 4\n'
    #com1 += 'LREAL = Auto; LPLANE = .TRUE.\n'
    #com1 += 'LSCALAPACK = .FALSE.\n\n'
    #f.write(com1)
    ###### 4: dft+D+U
    f.write('# functional (PE=PBE,RP=RPBE,RE=revPBE,b3=B3LYP,ML=vdw-df2, MK=rev-vdW-DF2, \n')
    f.write('# D correction if not vdW-DF (0-no, 1-d2, 11-d3_zero(Grimme), 12-de_BJ, 2-ts)\n')
    comm = "GGA = " + gga + "\n"
    ### B3LYP
    if re.search('b3',dic['dft'],re.IGNORECASE):
        comm += 'LHFCALC = .TRUE.\n'
        comm += 'ALGO = D; TIME = 0.5\n'
        comm += 'AEXX = 0.2\n'
        comm += 'AGGAX = 0.72\n'
        comm += 'AGGAC = 0.81\n'
        comm += 'ALDAC = 0.19\n'
    ### PBE style
    elif re.search('hs',dic['dft'],re.IGNORECASE):
        comm += 'LHFCALC = .TRUE.\n'
        comm += 'ALGO = D; TIME = 0.4\n'
        comm += 'HFSCREEN = 0.207\n'
        comm += 'PRECFOCK = F\n'
        comm += '#block the NPAR\n'
    elif re.search('pe',dic['dft'],re.IGNORECASE) or re.search('re',dic['dft'],re.IGNORECASE):      # for PBE, revPBE
        L_paw = 1
        comm += 'LREAL = Auto; LPLANE = .TRUE.\n'
        comm += 'LSCALAPACK = .FALSE.\n\n'
        if not L_hybrid:
            comm += 'ALGO = Fast\n'
            comm += '#IALGO=48\n'
            comm += 'NSIM = 4\n'
            comm += 'NPAR = 4\n'
        else:                                       # PBE0 or revPBE0
            comm += 'LHFCALC = .TRUE.\n'
            comm += 'ALGO = D; TIME = 0.4 # IALGO=53\n'
            comm += '#NSIM = 4\n'
            comm += 'ENCUTFOCK = 0\n'
            comm += 'NPAR = 4 ! NCORE*NPAR=total core; NCORE handles one band\n'
    elif re.search('e0',dic['dft'],re.IGNORECASE):
        comm += 'LHFCALC = .TRUE.\n'
        comm += 'ALGO = D; TIME = 0.4\n'
        comm += 'PRECFOCK = F\n'
        comm += '#block the NPAR\n'
    ### dispersion included dft?
    #if re.search('ML',dic['dft'],re.IGNORECASE) or re.search('MK',dic['dft'],re.IGNORECASE):
    if L_vdw:
        comm += "#vdW-DF parameter defined here: revPBE-DF, \n"
        comm += 'LUSE_VDW = .TRUE. \n'
        #comm += 'Zab_vdW = 1.8867 \n'
        comm += 'AGGAC = 0.0000 \n'
        #comm += 'LASPH = .TRUE. \n'
        if re.search('MK',dic['dft'],re.IGNORECASE):
            comm += 'PARAM1 = 0.1234 \n'
            comm += 'PARAM2 = 0.711357 \n'
    else:            
        if dic['dispersion']=='d2':
            comm += 'IVDW = 10      ! D2\n\n'
        elif dic['dispersion']=='d3':
            comm += 'IVDW = 11      ! D3-Grimme\n\n'
        elif dic['dispersion']=='d3bj':
            comm += 'IVDW = 12      ! D3-Becke-Jonson\n\n'
        elif dic['dispersion']=='ts':
            comm += 'IVDW = 20      ! Tkatchenko-Scheffler\n\n'
    comm += '\n'
    f.write(comm)
    ### DFT virtual orbital && read CHGCAR
    comm = '### DFT virtual orbital && CHGCAR\n'
    comm += '#ALGO = Exact\n'
    comm += '#NBANDS = 64\n'
    comm += '#LOPTICS = .TRUE.\n'
    if L_paw and re.match('1',str(icharg)):
        comm += 'LMAXMIX = 4 twice l of PP, s,p:2 l:4 f:6\n\n'
    else:
        comm += '#LMAXMIX = 4 for PAW && ICHARG = 1\n\n'
    f.write(comm)

    f.write('# GGA more\n')
    com1 = '#GGA_COMPAT=.FALSE. fpr bulk symmetry\n'
    com1 += '#VOSKOWN = 1 for PW91\n'
    f.write(com1)
    f.write('### U-correction\n')
    if dic['uterm']:
        ldaj = 1.0
        ldau = ldaj + dic['uterm']
        comm = 'LDAU = .TRUE.\nLDAUTYPE = 2\nLDAUL = 2 -1 -1 -1\nLDAUU = '+str(ldau)+' 0 0 0\nLDAJ = 1 0 0 0\n\n'
    else:
        comm = '\n'
        f.write(comm)
    ### 5: Movement: Relaxation, MD
    f.write('# Optimization\n')
    if dic['relax'] == 'sp':
        comm = '#NSW = ; ISIF = ; IBRION = ; POTIM = \n\n'
    else:
        if dic['relax'] == 'atom' :
            isif = 2
        else:
            isif = 3
        if L_md:
            isif=0
            nsw=1000
            ibrion=0
            potim=1.0
        else:
            nsw=999
            ibrion=2
            potim=0.3
        comm = 'NSW = %d ; ISIF = %d\n' % (nsw, isif)
        comm += 'IBRION = %d\nPOTIM = %.1f in femto second\n' % (ibrion, potim)
        comm += '### AIMD more\n'
        if L_md:                #dic['dynamics']:
            comm += 'TEBEG=%d; TEEND=%d\n' % (tebeg, teend)
            if tebeg == teend:
                smass = -3
            else:
                smass = -1
            comm += 'ISYM=0; SMASS = %d     for standard NVE ensemble w. const T \n\n' % smass
        else:
            if 'isym' in locals():
                comm += 'ISYM = 0'       
            comm += '#TEBEG = ; TEEND = \n' 
            comm += '#ISYM=0; SMASS = 0.05 ! -3\-1 for NVE, >=0 for NVT \n\n'
            
    f.write(comm)
    ### AIMD
    #f.write('# aimd - nve:: nsw=5000; isif=0; ibrion=0;potim=2.0;tein=300;tebeg=300;teend=300;smass=0.05;isym=0\n\n')
    ### Solvent effect
    f.write('# Solvent effect::higher Ecut is required (LSOL=.TRUE. EB_K for dielectric) \n')
    comm=""
    if dic['solvent']:
        comm = 'LSOL = .TRUE.\n'
        if dic['dielectric']:
            comm += 'EB_K = ' + dic['dielectric'] + '\n'
        comm += '#TAU = 0\n'
        comm += '#LRHOB = .TRUE\n\n'
        f.write(comm)
    ###### 6: POST-SCF 
    f.write('###### POST SCF CALC :: SOC [DOS|PCHG|NEB]  #########\n')
    ### SOC
    f.write('### spin-orbit coupling   (LSORBIT=.TRUE.) \n')
    if dic['soc']:
        comm = 'LSORBIT = .TRUE.\n\n'
        f.write(comm)
    ### DOS
    if (dic['postscf'] == 'dos'):
        comm = '### DOS    ( EMIN, EMAX, NEDOS=6000) \n'
        comm += 'EMIN = -40\n'
        comm += 'EMAX = 20\n'
        comm += 'NEDOS = 6000\n'
    ### PCHG  
    elif (dic['postscf'] == 'pchg'):
        comm = '### Partial charge  (NBANDS, LPARD IBAND, KPUSE, LSEPB, LSEPK)\n'
        comm += 'NBANDS = 112\n'
        comm += 'LPARD = .TRUE.\n'
        comm += 'IBAND = 91 92 \n'
        comm += 'KPUSE = 1 2 3 4 \n'
        comm += 'LSEPB = .TRUE. \n'
        comm += 'LSEPK = .TRUE. \n'
    ### NEB
    elif (dic['postscf'] == 'neb'):
        f.write('### NEB \n')

    f.write(comm)

    ###### 7: LOGFILE 
    f.write('###### log\n')
    if 0 < ilog:
        lwave = '.TRUE.'
    if 1 < ilog:
        lcharg = '.TRUE.'
    if 2 < ilog:
        laechg = '.TRUE.'
    comm = 'LROBIT = 11\nLAECHG = %s\nLWAVE = %s\nLCHARG = %s\n\n' % (laechg, lwave, lcharg)
    f.write(comm)

    ###### Extras
    comm = '#IDIPOL = 3\n'
    f.write(comm)

    f.close()
    return 0

def main():
    parser = argparse.ArgumentParser(description='Test function ')
    parser.add_argument('-j', '--job', choices=['getmag','ak','Ucorr'], help='read POSCAR then get MAGMOM')
    parser.add_argument('-s', '--poscar', default='POSCAR', help='POSCAR to be read')
    parser.add_argument('-ind', '--index', type=int, help='return index of Ucorr list')
    parser.add_argument('-m', '--magmom',  nargs='*', help='magmom input as hash')
    args = parser.parse_args()

    if args.job == 'getmag':
        st = make_mag_4pos(args.poscar, magin=args.magmom)
        print(st)
    elif args.job == 'ak':
        _, atoms = get_atoms_4pos(args.poscar)
        print(atoms)
    elif args.job == 'Ucorr':
        _, atoms = get_atoms_4pos(args.poscar)
        ### to be called in bash
        ldaul, ldauu, ldauj = get_Ucorr(atoms)
        if args.index == 0:
            print(ldaul)
            return ldaul
        elif args.index == 1:
            print(ldauu)
            return ldauu
        elif args.index == 2:
            print(ldauj)
            return ldauj
        else:
            print (ldaul, ldauu, ldauj)
            return ldaul, ldauu, ldauj
    return 0


if __name__ == '__main__':
    print(f"in main: {__name__}")
    main()
