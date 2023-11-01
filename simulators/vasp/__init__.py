from ...atoms   import *
from ...ncio    import cleansymb, get_unique_symbs, convert_xyz2abc
from ...units   import R
from ...units.sci_func import fermi_dirac, derivative_fd
from glob import glob
import os, math, sys
import numpy as np
from ...aux.aux_parse import parse_line, parse_lines
 
### import io.read, io.write inside Vasp class

# VASP Simulation Object
# made by Noh           2021. 8.
# modified by J. Park   2021.10. class Vasp, PDOS[nonmag]
#

class Vasp(object):

    """
    Vasp(atoms)

    Class for management of VASP simulation
    class varible
        nitem
        cpfiles     for copy after running simulator
        mvfiles     for mv after running simulator
        checkfile   to check whether checkfile exists for running simulator
        optfile     opt file to be read for vib calculation
    Parameters
    ----------
    symbol  : AtomsSystem
        Class instance of AtomsSystem

    Optional parameters
    -------------------
    """
    nitem = 0
    cpfiles = ['POSCAR', 'CONTCAR', 'OSZICAR']
    mvfiles = ['OUTCAR', 'XDATCAR']
    checkfile = 'OUTCAR'
    optfile = 'CONTCAR'

    def __init__(self, atoms):
        self.__class__.nitem += 1
        self.atoms = atoms 
        self._params = {
            # 1. Name and basic options       
            'SYSTEM'    :     'vasp',       # text, system name
            'NPAR'      :          1,       # integer, number of bands
            'NCORE'     :          0,       # integer, number of cores
            'IBRION'    :          2,       # 2=CG/Default, 5=Hessian 
            'LWAVE'     :        'F',       # boolean, write WAVECAR
            'LCHARG'    :        'F',       # boolean, write CHGCAR
            'NSW'       :          0,       # integer, optimization step
            'PREC'      : 'Accurate',       # precision (Low, Normal, Accurate)
            'ALGO'      :     'FAST',       # Algorithm, GGA/LDA=Normal, Fast
            'ISTART'    :          0,       # 0:new, 1:WAVECAR, 2:samecutoff
            'ICHARG'    :          2,       # charge from 0:WAVECAR, 1:file, 2:atomic, 11:keep CHARGE
            'ISIF'      :          2,       # 2=Constant cell, 3=relax cell
            # 2. SCF/kgrid/functional parameters          
            'ENCUT'     :        400,       # float, plane wave basis energy cutoff
            'ISMEAR'    :          0,       # integer, smearing 0=Gauss, 1=Metal
            'SIGMA'     :       0.05,       # positive float
            'NSIM'      :          1,       # integer, bands optimized in RMM-DIIS
            'NELMIN'    :          4,       # integer, min SCF steps
            'NELM'      :        500,       # integer, max SCF steps
            'EDIFF'     :     0.0001,       # float, tolerance of ground state
            'KPOINTS'   :  [1, 1, 1],       # list, 3-vector
            'XC'        :      'GGA',       # GGA, LDA
            'XCAUTHOR'  :      'PE' ,       # PE=PBE, 91=PW91, RP=Revised PBE 
            # 3. Optional parameters 
            'POTIM'     :        0.3,       # displacement  
            'EDIFFG'    :      -0.05,       # float, stopping relaxation loop
            'IVDW'      :         12,       # 11: D3 zero damping 12: D3 BJ damping
            'LDIPOL'    :        'F',       # dipole correction
            'IDIPOL'    :          3,       # 1: x, 2: y, 3: z, 4: all
            'LPLANE'    :        'T',       # data distribution over Nodes
            'ADDGRID'   :        'T',       # add grid for charge augmentation
            'LREAL'     :     'Auto',       # for real space projection
            'ISYM'      :         -1,       # -1 = symmetry off completely
            'LASPH'     :        'T',       # non-spherical contribtuion
            'LMAXMIX'   :          4,       # Density Mixer handles quantumNumber upto (4: d-elements, 6: f-elements)
            'ISPIN'     :          1,       # 1 = Spin-restricted, 2 = spin-unrestricted
            'MAGMOM'    :         None,       # can be in or use default
            }

    def get_options(self):
        """
        print the list of available options and their default values
       
        Parameters
        ----------

        Optional parameters
        -------------------

        Example
        -------
        >>> sim.get_options()
        """
        return self._params.items()
    ### useless: deprecate
    def read_file(fname):
        lineinfo = []
        wordinfo = []
        with open(fname) as f:
            for i, l in enumerate(f):
                line = l
                word = line.split()
                lineinfo.append(line)
                wordinfo.append(word)

        return lineinfo, wordinfo
    
    def set_option(self, key, value):
        
        """
        change the options
        available key and default values
        --------------------------------
        Parameters
        ----------
        key: str
            option name
        value: (various)
            option name
        Optional parameters
        -------------------
        Example
        -------
        >>> sim.set_options('KPOINTS', [5, 5, 1])
        """
        if key not in self._params.keys():
            raise ValueError("Invalid option," + key)
        else:
            self._params[key] = value
    def set_options(self, **kw):
        for k, v in kw.items():
            key = k.upper()
            ### this can be change or add
            self._params[key] = v
        return 0

    def set_atoms(self, contcar):
        self.atoms = readAtomicStructure(contcar)

    ### might be redundant with vasp.write_poscar
    def write_POSCAR(self, file_name='POSCAR', mode='cartesian', fix=None):
        components = self.atoms.get_contents().items()
        message  = ' '
        for i in components:
            message = message + str(i[0]) + '   '
        cell1    = self.atoms.get_cell()[0]
        cell2    = self.atoms.get_cell()[1]
        cell3    = self.atoms.get_cell()[2]

        #------------- POSCAR --------------------
        fout = open(file_name, 'w')
        fout.write("%s\n" % message)
        fout.write("1.000  # fixed lattice parameter unit\n")
        fout.write("%15.9f %15.9f %15.9f\n" % tuple(cell1))
        fout.write("%15.9f %15.9f %15.9f\n" % tuple(cell2))
        fout.write("%15.9f %15.9f %15.9f\n" % tuple(cell3))
        atm_line = ''; len_line = ''
        lines = []
        for sym, num in components:     # sym: symbol
            self.atoms.select_elements(sym)
            #atoms1 = self.atoms.copy_atoms() # sorted order is lost in a.copy_atoms()
            atm_line = atm_line + sym      + '   '
            len_line = len_line + str(num) + '   ' 
            for atom in self.atoms:
                x = 0. ; y = 0.; z = 0.
                if mode == 'cartesian':
                   x, y, z = Vector(atom.get_position())
                elif mode == 'direct':
                   x, y, z = Vector(atom.get_position())
                   x = x/(cell1[0] + cell1[1] + cell1[2])
                   y = y/(cell2[0] + cell2[1] + cell2[2])
                   z = z/(cell3[0] + cell3[1] + cell3[2])
                sline = f"{x:15.9f} {y:15.9f} {z:15.9f}"
                lines.append(sline)
                #print(sline)
        atm_line += '\n'; len_line += '\n'
        fout.write(atm_line)
        fout.write(len_line)
        fout.write("Selective Dynamics # constraints enabled\n")

        if mode == "cartesian":
            fout.write("Cartesian \n")
        elif mode == "direct":
            fout.write("Direct \n")
        
        for i in range(len(lines)):
            idx = i+1
            if fix == None:
                fout.write(str(lines[i]) + "   T   T   T \n")
            elif fix is not None:
                if idx in fix:
                    fout.write(str(lines[i]) + "   F   F   F \n")
                else:
                    fout.write(str(lines[i]) + "   T   T   T \n")
        fout.close()
    
    

    def write_KPOINTS(self):
        #-------------KPOINTS-------------------        
        p = self._params
        KPOINTS = open('KPOINTS', 'w')
        KPOINTS.write("k-points\n")
        KPOINTS.write("0\n")
        KPOINTS.write("G\n")
        KPOINTS.write("%i %i %i\n" % (p['KPOINTS'][0], p['KPOINTS'][1], p['KPOINTS'][2]))
        KPOINTS.write("0 0 0 \n")
        KPOINTS.close()

    def write_POTCAR(self, xc='PBE'):
        #-------------POTCAR--------------------
        p = self._params
        from ...aux.env import vasp_POTCAR_LDA  as LDA_path
        from ...aux.env import vasp_POTCAR_PBE  as PBE_path
        from ...aux.env import vasp_POTCAR_PW91 as PW91_path

        if xc == 'PBE':
            POTCAR_PATH = PBE_path
        elif xc == 'LDA':
            POTCAR_PATH = LDA_path
        elif xc == 'PW91':
            POTCAR_PATH = PW91_path
        else:
            print("select type of XC in PBE, LDA, PW91")

        components = self.atoms.get_contents().items()
        element = []; element_refine = []
        for sym, num in components:
            element.append(sym)
        
        sv_list = ["Li", "K", "Rb", "Cs", "Sr", "Ba", "Sc", "Y", "Zr"]
        pv_list = ["Na", "Ca", "Ti", "V", "Cr", "Mn", "Nb", "Mo", "Tc", "Hf", "Ta", "W", "Os"]
        d_list  = ["Ga", "Ge", "In", "Sn", "Tl", "Pb", "Bi", "Po", "At"]

        for name in element:
            if name in sv_list:
                element_refine.append(name+'_sv')
            elif name in pv_list:
                element_refine.append(name+'_pv')
            elif name in d_list:
                element_refine.append(name+'_d')
            else:
                element_refine.append(name)
        
        cmd = 'cat'
        os.system('rm -rf POTCAR')
        for element in element_refine:
            addcmd = ' ' + '%s/%s/POTCAR' % (POTCAR_PATH, element)
            cmd = cmd + addcmd
        cmd = cmd + ' > POTCAR'
        os.system('%s' % cmd)

    def write_INCAR(self):
        #-------------INCAR---------------------
        p = self._params
        #print(f"in writing INCAR {p['NPAR']}")
        INCAR = open('INCAR', 'w')
        INCAR.write("# VASP basic control parameters\n\n")
        INCAR.write("SYSTEM        =   %s\n" % p['SYSTEM'])
        INCAR.write("ISTART        =   %i\n" % p['ISTART']) 
        INCAR.write("ICHARG        =   %i\n" % p['ICHARG']) 
        INCAR.write("ISIF          =   %i\n\n" % p['ISIF']) 
        INCAR.write("IBRION        =   %i\n" % int(p['IBRION']))
        INCAR.write("NSW           =   %i\n" % p['NSW']) 
        INCAR.write("PREC          =   %s\n" % p['PREC']) 
        INCAR.write("ALGO          =   %s\n" % p['ALGO'])

        if p['NCORE'] :
            INCAR.write(f"{'NCORE':<15}={p['NCORE']:5d}\n")
        else:
            INCAR.write(f"{'NPAR':<15}={p['NPAR']:5d}\n")
        
        INCAR.write("LWAVE         =   %s\n" % p['LWAVE']) 
        INCAR.write("LCHARG        =   %s\n" % p['LCHARG']) 

        INCAR.write("# VASP convergence parameters \n\n")
        INCAR.write("ENCUT         =   %f\n" % float(p['ENCUT']))
        INCAR.write("ISMEAR        =   %i\n" % p['ISMEAR'])
        INCAR.write("SIGMA         =   %f\n" % p['SIGMA'])
        INCAR.write("NSIM          =   %i\n" % p['NSIM'])
        INCAR.write("NELMIN        =   %i\n" % p['NELMIN'])
        INCAR.write("NELM          =   %i\n" % p['NELM'])
        INCAR.write("EDIFF         =   %f\n" % float(p['EDIFF']))
        INCAR.write("EDIFFG        =   %f\n" % float(p['EDIFFG']))
        INCAR.write("%s           =   %s\n\n" % (p['XC'], p['XCAUTHOR']))
        INCAR.write("# VASP optional parameters \n\n")
        INCAR.write("POTIM         =   %f\n" % float(p['POTIM']))
        INCAR.write("IVDW          =   %i\n" % int(p['IVDW']))
        INCAR.write("LDIPOL        =   %s\n" % p['LDIPOL'])
        INCAR.write("IDIPOL        =   %i\n" % int(p['IDIPOL']))
        INCAR.write("LPLANE        =   %s\n" % p['LPLANE'])
        INCAR.write("ADDGRID       =   %s\n" % p['ADDGRID'])
        INCAR.write("LREAL         =   %s\n" % p['LREAL'])
        INCAR.write("ISYM          =   %i\n" % p['ISYM'])
        INCAR.write("LASPH         =   %s\n" % p['LASPH'])
        INCAR.write("LMAXMIX       =   %i\n" % p['LMAXMIX'])
        INCAR.write("ISPIN         =   %i\n\n" % p['ISPIN'])
        if p['ISPIN'] == 2:
            if p['MAGMOM']:
                #print(f"True {p['MAGMOM']}"), do no pass p['MAGMOM'] if it is {} -> it became True
                str_mag = get_magmom_4pos(pos="POSCAR", magin=f"{p['MAGMOM']}")
            else:
                str_mag = get_magmom_4pos(pos="POSCAR")
            INCAR.write(f"{str_mag}")
        INCAR.close()



    def run_catalysis(self, mode='sp', fix=None):
        """ 
        Run VASP with options
        mode    opt for optimization
                sp  
                vib for calc of vibrational frequency
        fix     atom index (starts from 1) to be fixed for vibration calc      
        """
        p = self._params
        
        #if p['SERVER'] == 'kisti':
        #    from oorinano.env.env_kisti import vasp_calculator as executable
        #else:
        from ...aux.env import vasp_calculator as executable


        ### obtain non-INCAR params
        if not 'KPOINTS' in p.keys():
            p['KPOINTS'] = [1, 1, 1]
        if 'NPROC' in p.keys():
            nproc = p['NPROC']
        else:
            nproc = 1
        del p['NPROC']

        ### these are set in set_params() 
        #p['NPAR']       = npar
        #p['ENCUT']      = encut
        #p['EDIFF']      = ediff
        #p['EDIFFG']     = ediffg
        #p['KPOINTS']    = kpoints
        
        if mode == 'opt':
            p['IBRION'] = 2
            p['POTIM']  = 0.300
            p['NSW']    = 500
        
        if mode == 'sp':
            p['IBRION'] = 2
            p['POTIM']  = 0.300
            p['NSW']    =   0

        if mode == 'vib':   # cal adsorbate for T*S
            p['IBRION'] = 5
            p['POTIM']  = 0.015
            p['NSW']    = 1

        # run_simulation
        cmd = f'mpirun -np {nproc}  {executable} > stdout.txt'

        self.write_POSCAR(fix=fix)
        self.write_KPOINTS()
        del p['KPOINTS']        # remove params not in INCAR
        self.write_POTCAR() 
        self.write_INCAR()      # remove p['SERVER']
        
        os.system(cmd)
    
    def save_files(self, fsuffix=None):
        for f in self.__class__.cpfiles:
            os.system(f'cp {f} {f}_{fsuffix}')
        for f in self.__class__.mvfiles:
            os.system(f'mv {f} {f}_{fsuffix}')
        return 0
    
    def save_checkfile(self, fsuffix=None):
        os.system(f'mv {self.__class__.checkfile} {self.__class__.checkfile}_{fsuffix}')

    def get_total_energy(self, output_name='OUTCAR'):
        
        line_info, word_info = Vasp.read_file(output_name)
        
        VASP_E = []
        for i in range(len(line_info)):
            if 'y  w' in line_info[i]:
                TE = float(word_info[i][-1])
                VASP_E.append(TE)
            else:
                pass
        
        min_E = min(VASP_E)

        return min_E 

    def get_vibration_energy(self, output_name='OUTCAR', Temp=298.15):
        """
        Example:
        --------
        from oorinano import vasp    
        ZPE, TS = vasp.get_vibration_energy(Temp=300)
        """
        
        line_info, word_info = Vasp.read_file(output_name)

        ZPE = 0.0; TS = 0.0
        RT  = R * Temp

        vib_vasp = []
        with open(output_name, "r") as fp:
            for line in fp:
                if 'THz' in line and not 'f/i' in line: # to remove imaginary freq
                    eles = line.split()
                    freq_E = 0.001*float(eles[-2])      # convert meV to eV
                    vib_vasp.append(freq_E)             # eV
                else:
                    pass
        
        for i in range(len(vib_vasp)):
            energy  = vib_vasp[i]
            x       = energy / RT
            v1      = x / (math.exp(x) - 1)
            vlog    = 1 - math.exp(-x)
            v2      = -math.log(vlog)

            #E_TS   = RT * (v1 + v2)
            E_TS   = RT * v2
            #print(f"Eentropy: freq {energy*1000:10.5f} : {RT*v1:10.5f} {RT*v2:10.5f}")
            ZPE = ZPE + 0.5*energy
            TS  = TS  + E_TS
        
        return ZPE, TS

    def get_vibration_spectrum(output_name='OUTCAR', start=0, end=6000, npts=None, width=20.0, matplot=1):               
        """
        Example:
        --------
        from oorinano import ncio       
        at = readAtomicStructure('POSCAR')
        at2 = vasp2.Vasp(at)
        at2.get_vibration_specctrum(output_name='OUTCAR_imag', matplot=1, start=-2000, end=6000)
        """
                                                                                                                          
        def read_file(fname):
            lineinfo = []
            wordinfo = []
            with open(fname) as f:
                for i, l in enumerate(f):
                    line = l
                    word = line.split()
                    lineinfo.append(line)
                    wordinfo.append(word)
         
            return lineinfo, wordinfo 
        
        line_info, word_info = read_file(output_name) 
                                                                                                                          
        vib_cm = []
        for i in range(len(line_info)):
            if 'THz' in line_info[i]:
                a = str(word_info[i][1].strip())
                if a[-1] == "=":
                    meV = float(word_info[i][-2])
                    convert = meV * -8.06554
                    vib_cm.append(convert)
                else:
                    meV = float(word_info[i][-2])
                    convert = meV * 8.06554
                    vib_cm.append(convert)
                                                                                                                          
        def fold(frequencies, intensities, start=0, end=6000, npts=None, width=20.0):
            if not npts:
                npts = int((end - start) / width * 10 + 1)
            prefactor = 1; sigma = width / 2. / np.sqrt(2. * np.log(2.))
            spectrum = np.empty(npts)
            energies = np.linspace(start, end, npts)
            for i, energy in enumerate(energies):
                energies[i] = energy
                spectrum[i] = (intensities * 0.5 * width / np.pi / ((frequencies - energy) ** 2 + 0.25 * width**2)).sum()
                                                                                                                          
            return [energies, prefactor * spectrum]
                                                                                                                          
        frequencies = vib_cm
        intensities = np.ones(len(frequencies))
        energies, spectrum = fold(frequencies, intensities, start=start, end=end, width=50.0)
        outdata = np.empty([len(energies), 2])
        outdata.T[0] = energies
        outdata.T[1] = spectrum
                                                                                                                          
        VDOS = open('VDOS.dat', 'w')
        for row in outdata:
            VDOS.write('%.3f %15.5e\n' % (row[0], row[1]))
        VDOS.close()

        if matplot:
            import matplotlib.pyplot as plt
            plt.figure(figsize=(7,5))
            plt.rcParams['axes.linewidth'] = 2
                                                                                
            plt.xticks(fontsize=10)
            plt.yticks(fontsize=10)
                                                                                
            line = plt.plot(energies, spectrum, linewidth='2', color='k')
            plt.xlabel('Frequency [cm$^{-1}$]', fontsize=20)
            plt.ylabel('VDOS [1/cm$^{-1}$]', fontsize=20)
            plt.savefig('VDOS.png', format='png', dpi=600, bbox_inches='tight')
        else:
            pass

### functions inside module vasp
def readAtomicStructure(file_name):
    f = open(file_name)
    lines = f.readlines()

    # system info.
    line_title = lines[0]
    line_cell_unit = float(lines[1].split()[0])
    line_cell1 = lines[2].split()
    line_cell2 = lines[3].split()
    line_cell3 = lines[4].split()
    line_symb  = lines[5].split()
    line_numb  =lines[6].split()

    # number of atoms
    n_system = 0
    for n in line_numb:
        n_system += int(n)

    # symbol list
    list_symb = []; index = 0
    for symb in line_symb:
        list_symb += [symb]*int(line_numb[index])
        index += 1

    # cell info.
    cell1 = []; cell2 = []; cell3 = []
    for v1 in line_cell1:
        cell1.append(line_cell_unit*float(v1))   
    for v2 in line_cell2:
        cell2.append(line_cell_unit*float(v2))    
    for v3 in line_cell3:
        cell3.append(line_cell_unit*float(v3))   
    cell = [cell1, cell2, cell3]

    # Constraint
    line_atoms = ''

    if lines[7].lower()[:9] == 'selective':
        line_atoms = lines[8:]
    else:
        line_atoms = lines[7:]

    i = 0
    for line in line_atoms:
        atoms = []
        i += 1

        # Cartesian        
        if line.split()[0][0:1].lower() == 'c':
            line_coord = line_atoms[i:i+n_system]
            j = 0
            for coord in line_coord:
                x ,y ,z = coord.split()[0], coord.split()[1], coord.split()[2]
                x = float(x); y = float(y); z = float(z)
                symb = list_symb[j]
                atoms.append(Atom(symb,(x,y,z)))
                j += 1
                #print x,y,z 
            atoms_obj = AtomsSystem(atoms)
            atoms_obj.set_cell(cell)
            #name = 'POSCAR.xyz'
            #io.write_xyz(name, atoms_obj)
            return atoms_obj

        # Factional            
        elif line.split()[0][0:1].lower() == 'd':
            line_coord = line_atoms[i:i+n_system]
            j = 0
            for coord in line_coord:
                xf,yf,zf = coord.split()[0], coord.split()[1], coord.split()[2]
                xf = float(xf); yf = float(yf); zf = float(zf)
                new_coord = xf*Vector(cell[0])+\
                            yf*Vector(cell[1])+\
                            zf*Vector(cell[2])
                x,y,z = new_coord[0], new_coord[1], new_coord[2]
                symb = list_symb[j]
                atoms.append(Atom(symb,(x,y,z)))
                j += 1
            atoms_obj = AtomsSystem(atoms)
            atoms_obj.set_cell(cell)
            #name = 'POSCAR.xyz'
            #io.write_xyz(name, atoms_obj)
            return atoms_obj


def writeAtomicStructure(atoms, file_name='POSCAR_xxyz', mode='cartesian', constraint=None):
    POSCAR = open(file_name, 'w')
    #POSCAR.write('%s\n' % params['title'])
    POSCAR.write('%s\n' % 'title')
    if atoms.get_cell() is not None:
        va,vb,vc = atoms.get_cell()
        va = Vector(va); vb = Vector(vb); vc = Vector(vc)
    else:
        raise ValueError("Cell info. is necessary.")
    POSCAR.write('1.000 # fixed lattice parameter unit\n')
    POSCAR.write("%15.9f %15.9f %15.9f\n" % tuple(va))
    POSCAR.write("%15.9f %15.9f %15.9f\n" % tuple(vb))
    POSCAR.write("%15.9f %15.9f %15.9f\n" % tuple(vc))
    contents = atoms.get_contents()
    print (contents)
    atm_line = ''; len_line = ''
    lines = []
    for sym, num in contents.items():
        atoms.select_elements(sym)
        atoms1 = atoms.copy_atoms()
        atm_line = atm_line + sym      + '  '
        len_line = len_line + str(num) + '  '
        for atom in atoms1:
            x = 0,; y = 0.; z = 0.
            if mode == 'cartesian':
                x,y,z = Vector(atom.get_position())
            elif mode == 'direct':
                x,y,z = Vector(atom.get_position())
                x= x/(va[0]+va[1]+va[2])
                y= y/(vb[0]+vb[1]+vb[2])
                z= z/(vc[0]+vc[1]+vc[2])
            #constraints?
            
            lines.append("%15.9f %15.9f %15.9f T T T\n" % (x,y,z))
            
            #lines.append("%15.9f %15.9f %15.9f\n" % (x,y,z))
    atm_line += '\n'; len_line += '\n'
    POSCAR.write(atm_line)
    POSCAR.write(len_line)
    POSCAR.write("Selective Dynamics # constraints enabled\n")

    if mode == 'cartesian':
        POSCAR.write("Cartesian # direct lattice\n")
    elif mode == 'direct':
        POSCAR.write("Direct # direct lattice\n")

    for line in lines:
        POSCAR.write(line)
    POSCAR.close()


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

def get_magmom_4pos(pos='POSCAR', magin=None):
    '''
	Used in making INCAR: 
    	Read POSCAR [and input magmom]
    	Return string of MAGMOM
    pos     POSCAR
    magin   MAGMOM input list of [ 'atom symbol', magmom, ... ]

    '''
    magmom = {}
    #print(f"{magin}")
    if magin:
        #print(f"{magin} is True")
        if type(magin) == list:
            li = iter(magin)
            magmom = dict(zip(li, li))
        elif type(magin) == dict:
            magmom = magin
    natoms, atoms = get_atoms_4pos(pos)
    magstr="MAGMOM = "
    Lmag = False
    for index, atom in enumerate(atoms):
        if magin and atom in magmom.keys():
            magstr += f"{natoms[index]}*{magmom[atom]*1.5} "
            Lmag = True
        elif atom in atom_prop.keys():
            magstr += f"{natoms[index]}*{atom_prop[atom][1]*1.5} "
            Lmag = True
        else:
            print("ERROR: no magmom in input and repo")
            sys.exit(10)
            #magstr += f"{natoms[index]}*0 "
    if Lmag: return magstr + "100*0"
    else:    return "# " + magstr

def calc_pdos(fname = 'DOSCAR', atom_list=[1], option=None):
    '''
    DOSCAR: 5 [preline] + (ngrid + 1 [energy headline]) * (natom + 1 [total DOS])
        ngrid   number of energy grid
    
    atom_list should starts from 1
    option: in case of spin = 2
        None    sum spin-up and spin-down
        split   write spin-up and spin-down with two columns
    '''

    nline_pre = 5
    with open(fname, 'r') as f:
        lines = f.readlines()
    #line_info, word_info = Vasp.read_file('DOSCAR')

    pre2d = parse_lines(lines[:5])

    ### prelude analysis
    natom = int(pre2d[0][0])

    ### parse energy title
    elist   = parse_line(lines[5])
    Emax    = float(elist[0])
    Emin    = float(elist[1])
    ngrid   = int(elist[2])                     # number of Ene grid
    E_fermi  = float(elist[3])

    ### test parsing of one line in energy block
    if len(parse_line(lines[6])) == 5:
        spin = 2    # spin polarized
    else:           # ncol = 3
        spin = 1
    
    ### Obtain TDOS = [nline_pre+1: nline_pre + ngrid +1] 

    Ei_f0   = []         # E w.r.t. E_fermi
    tdos    = []        # non-spin or spin-up
    tdos_cum = []
    if spin == 2:
        tdos_dn     = []
        tdos_cum_dn = []

    ### for TDOS
    for i, line in enumerate(lines[nline_pre+1:nline_pre+1+ngrid]):
        lele = parse_line(line)
        Ei_f0.append(float(lele[0])-E_fermi)
        tdos.append(float(lele[1]))
        if spin == 1:
            tdos_cum.append(float(lele[2]))       # Tdos total
        ### wrong order for spin==2, by J. Park
        else: # spin == 2:
            tdos_dn.append(float(lele[2]))    # this is TDOS-up-acc if spin == 2
            tdos_cum.append(float(lele[3]))
            tdos_cum_dn.append(float(lele[4]))    # TDOS-down-acc   if spin == 2
    
    with open('TDOS.dat', 'w') as f:
        for i in range(len(Ei_f0)):
            st = f"{Ei_f0[i]:11.3f}  {tdos[i]:10.4g}"
            if spin == 1:
                st += f"  {tdos_cum[i]:10.4g}"
            else: # spin == 2:
                st += f"  {tdos_dn[i]:10.4g}  {tdos_cum[i]:10.4g}  {tdos_cum_dn[i]:10.4g}"       # in the order of original DOSCAR
            f.write(st+"\n")
  

    ### 2D list [atom][energy] ###   
    atoms_s_up     = []   ;    atoms_s_dn     = []     
    atoms_py_up    = []   ;    atoms_py_dn    = []
    atoms_pz_up    = []   ;    atoms_pz_dn    = []
    atoms_px_up    = []   ;    atoms_px_dn    = []
    atoms_dxy_up   = []   ;    atoms_dxy_dn   = []
    atoms_dyz_up   = []   ;    atoms_dyz_dn   = []
    atoms_dz2_up   = []   ;    atoms_dz2_dn   = []
    atoms_dxz_up   = []   ;    atoms_dxz_dn   = []
    atoms_dx2_up   = []   ;    atoms_dx2_dn   = []
    
        
    ### scan atom_list:: Make 2D atoms list
    #print(f"atom list {atom_list}")
    for iatom in atom_list:
        ### Initialize for atom ###
        s_up     = []   ;    s_dn   = []
        py_up    = []   ;    py_dn  = []
        pz_up    = []   ;    pz_dn  = []
        px_up    = []   ;    px_dn  = []
        dxy_up   = []   ;    dxy_dn = []
        dyz_up   = []   ;    dyz_dn = []
        dz2_up   = []   ;    dz2_dn = []
        dxz_up   = []   ;    dxz_dn = []
        dx2_up   = []   ;    dx2_dn = []
        ### obtain block
        ipivot = nline_pre + (ngrid + 1) * (iatom)   # +1 (TDOS) -1 (at start)

        for i, line in enumerate(lines[ ipivot+1 : ipivot + 1 + ngrid]):
            item = parse_line(line)
            ### spin case
            if spin == 1: # only for the s,p,d + spin case
                s_up.append(float(item[1]))
                py_up.append(float(item[2])) ;  pz_up.append(float(item[3])) ;  px_up.append(float(item[4]))
                dxy_up.append(float(item[5]));  dyz_up.append(float(item[6]));  dz2_up.append(float(item[7]))
                dxz_up.append(float(item[8]));  dx2_up.append(float(item[9]))
            else:   # for spin == 2
                s_up.append(float(item[1]))   ; s_dn.append(float(item[2]))
                py_up.append(float(item[3]))  ; py_dn.append(float(item[4]));
                pz_up.append(float(item[5]))  ; pz_dn.append(float(item[6]))
                px_up.append(float(item[7]))  ; px_dn.append(float(item[8]))
                dxy_up.append(float(item[9])) ; dxy_dn.append(float(item[10]))
                dyz_up.append(float(item[11])); dyz_dn.append(float(item[12]))
                dz2_up.append(float(item[13])); dz2_dn.append(float(item[14]))
                dxz_up.append(float(item[15])); dxz_dn.append(float(item[16]))
                dx2_up.append(float(item[17])); dx2_dn.append(float(item[18]))
        ### one atom dos in the list of energy
        atoms_s_up.append(s_up); 
        atoms_py_up.append(py_up)   ;   atoms_pz_up.append(pz_up)   ; atoms_px_up.append(px_up)     
        atoms_dxy_up.append(dxy_up) ;   atoms_dyz_up.append(dyz_up) ; atoms_dz2_up.append(dz2_up)      
        atoms_dxz_up.append(dxz_up) ;   atoms_dx2_up.append(dx2_up)
        if spin == 2:    
            atoms_s_dn.append(s_dn)    
            atoms_py_dn.append(py_dn)   ;   atoms_pz_dn.append(pz_dn)   ; atoms_px_dn.append(px_dn)  
            atoms_dxy_dn.append(dxy_dn) ;   atoms_dyz_dn.append(dyz_dn) ; atoms_dz2_dn.append(dz2_dn) 
            atoms_dxz_dn.append(dxz_dn) ;   atoms_dx2_dn.append(dx2_dn)
    ### you might extract dos for each atom

    ### summation for atoms :: Make 1D lsum list
    lsum_s_up   = []   ;   lsum_s_dn   = []
    lsum_py_up  = []   ;   lsum_py_dn  = []
    lsum_pz_up  = []   ;   lsum_pz_dn  = []
    lsum_px_up  = []   ;   lsum_px_dn  = []
    lsum_dxy_up = []   ;   lsum_dxy_dn = []
    lsum_dyz_up = []   ;   lsum_dyz_dn = []
    lsum_dz2_up = []   ;   lsum_dz2_dn = []
    lsum_dxz_up = []   ;   lsum_dxz_dn = []
    lsum_dx2_up = []   ;   lsum_dx2_dn = []
    ### shape of atoms_dos=[natoms_in_list][ngrid]
    for j in range(ngrid):
        sum_s_up   = 0   ;   sum_s_dn   = 0
        sum_py_up  = 0   ;   sum_py_dn  = 0
        sum_pz_up  = 0   ;   sum_pz_dn  = 0
        sum_px_up  = 0   ;   sum_px_dn  = 0
        sum_dxy_up = 0   ;   sum_dxy_dn = 0
        sum_dyz_up = 0   ;   sum_dyz_dn = 0
        sum_dz2_up = 0   ;   sum_dz2_dn = 0
        sum_dxz_up = 0   ;   sum_dxz_dn = 0
        sum_dx2_up = 0   ;   sum_dx2_dn = 0
        for i in range(len(atom_list)):
            # atom_list index should starts from 1, len(s_up) == nEne * 'sum of Atoms'
            sum_s_up    += atoms_s_up[i][j]       
            sum_py_up   += atoms_py_up[i][j]     
            sum_pz_up   += atoms_pz_up[i][j]     
            sum_px_up   += atoms_px_up[i][j]     
            sum_dxy_up  += atoms_dxy_up[i][j]   
            sum_dyz_up  += atoms_dyz_up[i][j]   
            sum_dz2_up  += atoms_dz2_up[i][j]   
            sum_dxz_up  += atoms_dxz_up[i][j]   
            sum_dx2_up  += atoms_dx2_up[i][j]  
            if spin == 2:
                sum_s_dn    += atoms_s_dn[i][j]   
                sum_py_dn   += atoms_py_dn[i][j]  
                sum_pz_dn   += atoms_pz_dn[i][j]  
                sum_px_dn   += atoms_px_dn[i][j]  
                sum_dxy_dn  += atoms_dxy_dn[i][j] 
                sum_dyz_dn  += atoms_dyz_dn[i][j] 
                sum_dz2_dn  += atoms_dz2_dn[i][j] 
                sum_dxz_dn  += atoms_dxz_dn[i][j]                  
                sum_dx2_dn  += atoms_dx2_dn[i][j]              # each summed orbitals
        lsum_s_up.append(sum_s_up)   
        lsum_py_up.append(sum_py_up)  ; lsum_pz_up.append(sum_pz_up)  ; lsum_px_up.append(sum_px_up)
        lsum_dxy_up.append(sum_dxy_up); lsum_dyz_up.append(sum_dyz_up); lsum_dz2_up.append(sum_dz2_up)
        lsum_dxz_up.append(sum_dxz_up); lsum_dx2_up.append(sum_dx2_up)
        if spin == 2:
            lsum_s_dn.append(sum_s_dn)   
            lsum_py_dn.append(sum_py_dn)  ; lsum_pz_dn.append(sum_pz_dn)  ; lsum_px_dn.append(sum_px_dn)
            lsum_dxy_dn.append(sum_dxy_dn); lsum_dyz_dn.append(sum_dyz_dn); lsum_dz2_dn.append(sum_dz2_dn)
            lsum_dxz_dn.append(sum_dxz_dn); lsum_dx2_dn.append(sum_dx2_dn)

    ### make string for atom index
    if len(atom_list) == 1:
        froot = f"Atoms{atom_list[0]}"
    else:
        froot = f"Atoms{atom_list[0]}-{atom_list[-1]}_N{len(atom_list)}"
    if option:
        froot += option[:3]
    outfile = froot + '.dat'
    outf = open(outfile, 'w')
    lsum_p = []
    lsum_d = []
    lsum_p_dn = []
    lsum_d_dn = []
    #print(f"dim ngrid {range(ngrid)}, {len(lsum_px_up)}")
    for i in range(ngrid):
        lsum_p.append(lsum_px_up[i]+lsum_py_up[i]+lsum_pz_up[i])
        lsum_d.append(lsum_dx2_up[i]+lsum_dxy_up[i]+lsum_dxz_up[i]+lsum_dyz_up[i]+lsum_dz2_up[i])
        if spin == 1:
            ### format:: E, lsum_s, lsum_p, lsum_d, lsum(s+p+d)
            st = "%15.8f %15.8f %15.8f %15.8f %15.8f" % (Ei_f0[i],lsum_s_up[i], lsum_p[i], lsum_d[i], lsum_s_up[i]+lsum_p[i]+lsum_d[i])
        else:
            lsum_p_dn.append(lsum_px_dn[i]+lsum_py_dn[i]+lsum_pz_dn[i])
            lsum_d_dn.append(lsum_dx2_dn[i]+lsum_dxy_dn[i]+lsum_dxz_dn[i]+lsum_dyz_dn[i]+lsum_dz2_dn[i])
            if option == 'split':
                ### format:: E, lsum_s_up, lsum_s_dn, lsum_p, lsum_p_dn, lsum_d, lsum_d_dn, lsum(s+p+d), lsum(s+p+d)_dn 
                st = "%15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f" % \
                    (Ei_f0[i],lsum_s_up[i],lsum_s_dn[i],lsum_p[i],lsum_p_dn[i],lsum_d[i],lsum_d_dn[i],lsum_s_up[i]+lsum_p[i]+lsum_d[i],lsum_s_dn[i]+lsum_p_dn[i]+lsum_d_dn[i])
            elif option == 'polar':
                st = "%15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f" % \
                    (Ei_f0[i],lsum_s_up[i],-lsum_s_dn[i],lsum_p[i],-lsum_p_dn[i],lsum_d[i],-lsum_d_dn[i],lsum_s_up[i]+lsum_p[i]+lsum_d[i],-(lsum_s_dn[i]+lsum_p_dn[i]+lsum_d_dn[i]))           
            ### No option: lsum up and down ### format:: E, lsum_s, lsum_p, lsum_d, lsum(s+p+d)
            else:
                st = "%15.8f %15.8f %15.8f %15.8f %15.8f" % (Ei_f0[i],lsum_s_up[i]+lsum_s_dn[i],lsum_p[i]+lsum_p_dn[i],lsum_d[i]+lsum_d_dn[i],lsum_s_up[i]+lsum_p[i]+lsum_d[i]+lsum_s_dn[i]+lsum_p_dn[i]+lsum_d_dn[i])
        
        outf.write(st)
        outf.write("\n")
    outf.close()
    os.system(f'cp {fname} SUM_ATOM.dat')


def pdos_split_sum(atom_list=[1]):
    '''
    atom_list should starts from 1
    '''
    nline_pre = 6
    line_info, word_info = Vasp.read_file('DOSCAR')
    nEne = int(word_info[5][2])                     # number of E points
    E_fermi  = float(word_info[5][3])
    if len(word_info[6]) == 5:
        spin = 2    # spin polarized
    else:
        spin = 1
    #print(f"spin: {spin}, 2 for spin-polarized")

    ### Obtain TDOS

    E_f0   = []         # E w.r.t. E_fermi
    tdos    = []        # non-spin or spin-up
    tdos_cum = []
    if spin == 2:
        tdos_dn     = []
        tdos_cum_dn = []

    ### for TDOS
    for ind in range(nline_pre, nline_pre+nEne+1):  # ind: index of line in DOSCAR
        E_f0.append(float(word_info[ind][0])-E_fermi)
        tdos.append(float(word_info[ind][1]))
        if spin != 2:
            tdos_cum.append(float(word_info[ind][2]))       # this is TDOS-down   if spin == 2
        ### wrong order for spin==2, by J. Park
        else: # spin == 2:
            tdos_dn.append(float(word_info[ind][2]))    # this is TDOS-up-acc if spin == 2
            tdos_cum.append(float(word_info[ind][3]))
            tdos_cum_dn.append(float(word_info[ind][4]))    # TDOS-down-acc   if spin == 2

    
    with open('TDOS.dat', 'w') as f:
        for i in range(len(E_f0)):
            st = f"{E_f0[i]:11.3f}  {tdos[i]:10.4g}"
            if spin != 2:
                st += f"  {tdos_cum[i]:10.4g}"
            else: # spin == 2:
                st += f"  {tdos_dn[i]:10.4g}  {tdos_cum[i]:10.4g}  {tdos_cum_dn[i]:10.4g}"       # in the order of original DOSCAR
            f.write(st+"\n")


    ### Data list (s, p, d + spin) ###
    s_up     = []   ;    s_dn   = []
    py_up    = []   ;    py_dn  = []
    pz_up    = []   ;    pz_dn  = []
    px_up    = []   ;    px_dn  = []
    dxy_up   = []   ;    dxy_dn = []
    dyz_up   = []   ;    dyz_dn = []
    dz2_up   = []   ;    dz2_dn = []
    dxz_up   = []   ;    dxz_dn = []
    dx2_up   = []   ;    dx2_dn = []
    s_up_sum = []   ;  s_dn_sum = []
    p_up_sum = []   ;  p_dn_sum = []
    d_up_sum = []   ;  d_dn_sum = []
    ###

    E_f0i   = []
    
    ### scan all the line in DOSCAR
    for i in range(len(word_info)):
        ### spin case
        if len(word_info[i]) == 19: # only for the s,p,d + spin case
            item = word_info[i]
            E_f0i.append(float(item[0])-E_fermi)
            s_up.append(float(item[1]))
            s_dn.append(-float(item[2]))
            py_up.append(float(item[3]))
            py_dn.append(-float(item[4]))
            pz_up.append(float(item[5]))
            pz_dn.append(-float(item[6]))
            px_up.append(float(item[7]))
            px_dn.append(-float(item[8]))
            dxy_up.append(float(item[9]))
            dxy_dn.append(-float(item[10]))
            dyz_up.append(float(item[11]))
            dyz_dn.append(-float(item[12]))
            dz2_up.append(float(item[13]))
            dz2_dn.append(-float(item[14]))
            dxz_up.append(float(item[15]))
            dxz_dn.append(-float(item[16]))
            dx2_up.append(float(item[17]))
            dx2_dn.append(-float(item[18]))
        elif len(word_info[i]) == 10:
            item = word_info[i]
            E_f0i.append(float(item[0])-E_fermi)
            s_up.append(float(item[1]))
            py_up.append(float(item[2]))
            pz_up.append(float(item[3]))
            px_up.append(float(item[4]))
            dxy_up.append(float(item[5]))
            dyz_up.append(float(item[6]))
            dz2_up.append(float(item[7]))
            dxz_up.append(float(item[8]))
            dx2_up.append(float(item[9]))

    N_atoms = len(E_f0i) / nEne # double check
    
    
    ### in all atom split, index starts from 0
    for i in range(int(N_atoms)):
        dataDOS = open('%03d_ATOM.dat' % int(i+1), 'w')
        for j in range(nEne):
            if spin == 2:
                ### format:: E, sum_up, sum_dn, s_up, s_down, psum_up, psum_dn, dsum_up, dsum_dn, 18-lm: total 27 col
                st = "%15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f" % (E_f0i[i*nEne+j],s_up[i*nEne+j]+py_up[i*nEne+j]+pz_up[i*nEne+j]+px_up[i*nEne+j]+dxy_up[i*nEne+j]+dyz_up[i*nEne+j]+dz2_up[i*nEne+j]+dxz_up[i*nEne+j]+dx2_up[i*nEne+j],s_dn[i*nEne+j]+py_dn[i*nEne+j]+pz_dn[i*nEne+j]+px_dn[i*nEne+j]+dxy_dn[i*nEne+j]+dyz_dn[i*nEne+j]+dz2_dn[i*nEne+j]+dxz_dn[i*nEne+j]+dx2_dn[i*nEne+j],s_up[i*nEne+j],s_dn[i*nEne+j],py_up[i*nEne+j]+pz_up[i*nEne+j]+px_up[i*nEne+j],py_dn[i*nEne+j]+pz_dn[i*nEne+j]+px_dn[i*nEne+j],dxy_up[i*nEne+j]+dyz_up[i*nEne+j]+dz2_up[i*nEne+j]+dxz_up[i*nEne+j]+dx2_up[i*nEne+j],dxy_dn[i*nEne+j]+dyz_dn[i*nEne+j]+dz2_dn[i*nEne+j]+dxz_dn[i*nEne+j]+dx2_dn[i*nEne+j],s_up[i*nEne+j],s_dn[i*nEne+j],py_up[i*nEne+j],py_dn[i*nEne+j],pz_up[i*nEne+j],pz_dn[i*nEne+j],px_up[i*nEne+j],px_dn[i*nEne+j],dxy_up[i*nEne+j],dxy_dn[i*nEne+j],dyz_up[i*nEne+j],dyz_dn[i*nEne+j],dz2_up[i*nEne+j],dz2_dn[i*nEne+j],dxz_up[i*nEne+j],dxz_dn[i*nEne+j],dx2_up[i*nEne+j],dx2_dn[i*nEne+j])    
            else:
                st = "%15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f" % (E_f0i[i*nEne+j],s_up[i*nEne+j]+py_up[i*nEne+j]+pz_up[i*nEne+j]+px_up[i*nEne+j]+dxy_up[i*nEne+j]+dyz_up[i*nEne+j]+dz2_up[i*nEne+j]+dxz_up[i*nEne+j]+dx2_up[i*nEne+j],s_up[i*nEne+j],py_up[i*nEne+j]+pz_up[i*nEne+j]+px_up[i*nEne+j],dxy_up[i*nEne+j]+dyz_up[i*nEne+j]+dz2_up[i*nEne+j]+dxz_up[i*nEne+j]+dx2_up[i*nEne+j],s_up[i*nEne+j],py_up[i*nEne+j],pz_up[i*nEne+j],px_up[i*nEne+j],dxy_up[i*nEne+j],dyz_up[i*nEne+j],dz2_up[i*nEne+j],dxz_up[i*nEne+j],dx2_up[i*nEne+j])    
            dataDOS.write(st)
            dataDOS.write("\n")
        dataDOS.close()
    if not os.path.isdir('PDOS'):
        os.system('mkdir PDOS')
    os.system('mv *_ATOM.dat PDOS')
    
    ### Sum data list (s, p, d + spin) ###   
    sum_s_up     = []   ;    sum_s_dn     = []     
    sum_py_up    = []   ;    sum_py_dn    = []
    sum_pz_up    = []   ;    sum_pz_dn    = []
    sum_px_up    = []   ;    sum_px_dn    = []
    sum_dxy_up   = []   ;    sum_dxy_dn   = []
    sum_dyz_up   = []   ;    sum_dyz_dn   = []
    sum_dz2_up   = []   ;    sum_dz2_dn   = []
    sum_dxz_up   = []   ;    sum_dxz_dn   = []
    sum_dx2_up   = []   ;    sum_dx2_dn   = []
    sum_s_up_sum = []   ;    sum_s_dn_sum = []
    sum_p_up_sum = []   ;    sum_p_dn_sum = []
    sum_d_up_sum = []   ;    sum_d_dn_sum = []
    ###
    
    for j in range(nEne):
        num_s_up   = 0   ;   num_s_dn   = 0
        num_py_up  = 0   ;   num_py_dn  = 0
        num_pz_up  = 0   ;   num_pz_dn  = 0
        num_px_up  = 0   ;   num_px_dn  = 0
        num_dxy_up = 0   ;   num_dxy_dn = 0
        num_dyz_up = 0   ;   num_dyz_dn = 0
        num_dz2_up = 0   ;   num_dz2_dn = 0
        num_dxz_up = 0   ;   num_dxz_dn = 0
        num_dx2_up = 0   ;   num_dx2_dn = 0
        num_s_up_sum = 0 ;   num_s_dn_sum = 0
        num_p_up_sum = 0 ;   num_p_dn_sum = 0
        num_d_up_sum = 0 ;   num_d_dn_sum = 0
        for i in atom_list:
            # atom_list index should starts from 1, len(s_up) == nEne * 'num of Atoms'
            num_s_up += s_up[(i-1)*nEne+j]       
            num_py_up += py_up[(i-1)*nEne+j]     
            num_pz_up += pz_up[(i-1)*nEne+j]     
            num_px_up += px_up[(i-1)*nEne+j]     
            num_dxy_up += dxy_up[(i-1)*nEne+j]   
            num_dyz_up += dyz_up[(i-1)*nEne+j]   
            num_dz2_up += dz2_up[(i-1)*nEne+j]   
            num_dxz_up += dxz_up[(i-1)*nEne+j]   
            num_dx2_up += dx2_up[(i-1)*nEne+j]  
            if spin == 2:
                num_s_dn += s_dn[(i-1)*nEne+j]
                num_py_dn += py_dn[(i-1)*nEne+j]
                num_pz_dn += pz_dn[(i-1)*nEne+j]
                num_px_dn += px_dn[(i-1)*nEne+j]
                num_dxy_dn += dxy_dn[(i-1)*nEne+j]
                num_dyz_dn += dyz_dn[(i-1)*nEne+j]
                num_dz2_dn += dz2_dn[(i-1)*nEne+j]
                num_dxz_dn += dxz_dn[(i-1)*nEne+j]                 
                num_dx2_dn += dx2_dn[(i-1)*nEne+j]             # each summed orbitals
            num_s_up_sum += s_up[(i-1)*nEne+j]   
            num_p_up_sum += py_up[(i-1)*nEne+j] + pz_up[(i-1)*nEne+j] + px_up[(i-1)*nEne+j]
            num_d_up_sum += dxy_up[(i-1)*nEne+j] + dyz_up[(i-1)*nEne+j] + dz2_up[(i-1)*nEne+j] + dxz_up[(i-1)*nEne+j] + dx2_up[(i-1)*nEne+j]
            if spin == 2:
                num_s_dn_sum += s_dn[(i-1)*nEne+j]
                num_p_dn_sum += py_dn[(i-1)*nEne+j] + pz_dn[(i-1)*nEne+j] + px_dn[(i-1)*nEne+j]
                num_d_dn_sum += dxy_dn[(i-1)*nEne+j] + dyz_dn[(i-1)*nEne+j] + dz2_dn[(i-1)*nEne+j] + dxz_dn[(i-1)*nEne+j] + dx2_dn[(i-1)*nEne+j]
        # summation values: energy list
        sum_s_up.append(num_s_up)         ;   sum_s_dn.append(num_s_dn)
        sum_py_up.append(num_py_up)       ;   sum_py_dn.append(num_py_dn)
        sum_pz_up.append(num_pz_up)       ;   sum_pz_dn.append(num_pz_dn)
        sum_px_up.append(num_px_up)       ;   sum_px_dn.append(num_px_dn)
        sum_dxy_up.append(num_dxy_up)     ;   sum_dxy_dn.append(num_dxy_dn)
        sum_dyz_up.append(num_dyz_up)     ;   sum_dyz_dn.append(num_dyz_dn)
        sum_dz2_up.append(num_dz2_up)     ;   sum_dz2_dn.append(num_dz2_dn)
        sum_dxz_up.append(num_dxz_up)     ;   sum_dxz_dn.append(num_dxz_dn)
        sum_dx2_up.append(num_dx2_up)     ;   sum_dx2_dn.append(num_dx2_dn)
        sum_s_up_sum.append(num_s_up_sum) ;   sum_s_dn_sum.append(num_s_dn_sum)
        sum_p_up_sum.append(num_p_up_sum) ;   sum_p_dn_sum.append(num_p_dn_sum)
        sum_d_up_sum.append(num_d_up_sum) ;   sum_d_dn_sum.append(num_d_dn_sum)
    ### make string for atom index
    fname = f"Atoms{atom_list[0]}_{atom_list[-1]}_N{len(atom_list)}.dat"
    dataSUMDOS = open(fname, 'w')
    for i in range(nEne):
        if spin == 2:
            ### format:: E, sum_up, sum_dn, s_up, s_down, psum_up, psum_dn, dsum_up, dsum_dn, 18-lm: total 27 col
            st = "%15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f" % (E_f0i[i],sum_s_up_sum[i]+sum_p_up_sum[i]+sum_d_up_sum[i],sum_s_dn_sum[i]+sum_p_dn_sum[i]+sum_d_dn_sum[i],sum_s_up_sum[i],sum_s_dn_sum[i],sum_p_up_sum[i],sum_p_dn_sum[i],sum_d_up_sum[i],sum_d_dn_sum[i],sum_s_up[i],sum_s_dn[i],sum_py_up[i],sum_py_dn[i],sum_pz_up[i],sum_pz_dn[i],sum_px_up[i],sum_px_dn[i],sum_dxy_up[i],sum_dxy_dn[i],sum_dyz_up[i],sum_dyz_dn[i],sum_dz2_up[i],sum_dz2_dn[i],sum_dxz_up[i],sum_dxz_dn[i],sum_dx2_up[i],sum_dx2_dn[i])
        else:
            st = "%15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f" % (E_f0i[i],sum_s_up_sum[i]+sum_p_up_sum[i]+sum_d_up_sum[i],sum_s_up_sum[i],sum_p_up_sum[i],sum_d_up_sum[i],sum_s_up[i],sum_py_up[i],sum_pz_up[i],sum_px_up[i],sum_dxy_up[i],sum_dyz_up[i],sum_dz2_up[i],sum_dxz_up[i],sum_dx2_up[i])
        dataSUMDOS.write(st)
        dataSUMDOS.write("\n")
    dataSUMDOS.close()
    os.system(f'cp {fname} PDOS')
    os.system(f'cp {fname} SUM_ATOM.dat')

def pdos_orbital_analysis(fname='SUM_ATOM.dat', orbitals=1):
    '''
    d-band center theory : Nature, 376, 238 (1995)
    Fermi abundance      : J. Phys. Chem. C 121, 1530 (2017)
    highest peak         : Nat. Energy, 1, 16130 (2016)
    orbitals    0 for 's', 1 for 'p', 2 for 'd'
    '''
    import os, sys, math
    import numpy as np
    from scipy.misc import derivative
 
    def find_index(data, target):
        res = []
        lis = data
        while True:
            try:
                res.append(lis.index(target) + (res[-1]+1 if len(res)!=0 else 0))
                lis = data[res[-1]+1:]
            except:
                break
        return res
    
    ### Start of calculation
    line_info, word_info = Vasp.read_file(fname)
    
    E = [] ; Orb = []
    
    ### Obtain E and Orb from ortibals(s=0, p=1, etc)
    for i in range(len(word_info)):
        E.append(float(word_info[i][0]))
        Orb.append(float(word_info[i][orbitals]))

    for i in range(len(E)-1):
        if E[i] * E[i+1] <= 0:
            Fermi_idx = i + 1   # i < E_fermi < i+1
            break
    
    # find parameteres
    dE = E[1] - E[0]
    n_Orb = 0 ; E_n_Orb = 0
    df_n_Orb = 0 ; df_E_n_Orb = 0
    
    Orb_Ef = Orb[:Fermi_idx]
    if sum(Orb_Ef) >= 0:
        Orb_max = max(Orb_Ef) 
    elif sum(Orb_Ef) < 0:
        Orb_max = min(Orb_Ef)
    else:
        print("The data is not consistent for spin up or dn")
    
    E_Orb_max = find_index(Orb_Ef, Orb_max)
    print(f"density max indes {E_Orb_max}")

    if len(E_Orb_max) == 1:
        E_max = E[E_Orb_max[0]]
    else:
        print("Same max values exist")

    # calculation 
    for i in range(Fermi_idx):
        n_Orb       = n_Orb      + dE * Orb_Ef[i]
        E_n_Orb     = E_n_Orb    + dE * Orb_Ef[i] * E[i]
        df_n_Orb    = df_n_Orb   + dE * Orb_Ef[i] * derivative_fd(E[i])
        df_E_n_Orb  = df_E_n_Orb + dE * Orb_Ef[i] * E[i] * derivative_fd(E[i])

    Orb_cent = E_n_Orb / n_Orb
    Ef_abund = df_E_n_Orb / df_n_Orb
    print('orbital center ', '%12.8f' % Orb_cent)
    print('fermi-abudnace ', '%12.8f' % Ef_abund)
    print('highest peak at','%12.8f'  %  E_max, 'eV has', '%12.8f' % Orb_Ef[E_Orb_max[0]], 'states')
    
    return Orb_cent, Ef_abund, (E_max, Orb_Ef[E_Orb_max[0]])



