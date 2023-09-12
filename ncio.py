#
# read/write structure files
#

from __future__ import print_function
import re, sys
from math import sin, cos, sqrt, pi
from .units import ang2bohr, degrad
#from .atoms import atomic_weight, atomic_symbol, atomic_number
from .atoms import *

def write_struct(atoms, cellparameter=1.0, fname = "STRUCT.fdf"):

    if atoms.get_cell() != 'None':
        cell1 = atoms.get_cell()[0]
        cell2 = atoms.get_cell()[1]
        cell3 = atoms.get_cell()[2]

    #---------------STRUCT.fdf----------------
    fileS = open(fname, 'w')
    natm = len(atoms)
    fileS.write("NumberOfAtoms    %d           # Number of atoms\n" % natm)
    unique_symbs = get_unique_symbs(atoms)
    fileS.write("NumberOfSpecies  %d           # Number of species\n\n" % len(unique_symbs))
    fileS.write("%block ChemicalSpeciesLabel\n")

    for symb in unique_symbs:
        sym = ''.join(re.findall('[a-zA-Z]', symb))
        fileS.write(" %d %d %s\n" % (unique_symbs.index(symb)+1,atomic_number(sym),symb) )
    fileS.write("%endblock ChemicalSpeciesLabel\n")

    #Lattice
    fileS.write("\n#(3) Lattice, coordinates, k-sampling\n\n")
    fileS.write("LatticeConstant   %15.9f Ang\n" % cellparameter)
    fileS.write("%block LatticeVectors\n")
    if atoms.get_cell() != 'None':
        va, vb, vc = cell1, cell2, cell3
        fileS.write("%15.9f %15.9f %15.9f\n" % tuple(va))
        fileS.write("%15.9f %15.9f %15.9f\n" % tuple(vb))
        fileS.write("%15.9f %15.9f %15.9f\n" % tuple(vc))
    fileS.write("%endblock LatticeVectors\n\n")

    #Coordinates
    fileS.write("AtomicCoordinatesFormat Ang\n")
    fileS.write("%block AtomicCoordinatesAndAtomicSpecies\n")

    for atom in atoms:
        x,y,z = atom.get_position(); symb = atom.get_symbol()
        if atom.get_groupid() != 0:
            symb = symb + str(atom.get_groupid())
        fileS.write(" %15.9f %15.9f %15.9f %4d %4d\n" %\
                    (x,y,z,unique_symbs.index(symb)+1, atom.get_serial()))
        
    fileS.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")
    fileS.close()

def read_struct(file_name):
    vec_block = []; atoms_block = []; abc_cell_block = []
    atoms_length = 0; species = []
    n_of_species = 0; name = ''; atoms = []; cell = []; cell_scale = ''
    lattice_constant = 0.
    _is_ang_scale = 0; _is_bohr_scale = 0; _is_scaled_ang_scale = 0
    _is_fraction_scale = 0

    f = open(file_name)
    lines = f.readlines()

    i = 0
    for line in lines:
        #print i
        
        line_s = line.split(); keyword = ''

        if line_s:
            keyword = line_s[0].lower()
            #print keyword

        if keyword == "systemlabel":
            name = line_s[1]

        elif keyword == "latticeconstant":
            lattice_constant = float(line_s[1])
            #print lattice_constant
            try:
                cell_scale = line_s[2]
            except:
                cell_scale = 'Ang'

        elif keyword == "atomiccoordinatesformat":
            if line_s[1].lower() == 'ang':
                _is_ang_scale = 1
            elif line_s[1].lower() == 'bohr':
                _is_bohr_scale = 1
            elif line_s[1].lower() == 'scaledcartesian':
                #print "ON"
                _is_scaled_ang_scale = 1
            elif line_s[1].lower() == 'fractional':
                _is_fraction_scale = 1
            else:
                #print 'Warning : Default atomic scale, "Ang".\n'
                pass

        elif keyword == "numberofatoms":
            atoms_length = int(line_s[1])
            #print "natms", atoms_length

        elif keyword == "numberofspecies":
            n_of_species = int(line_s[1])
            #print "nspec", n_of_species

        elif keyword =="%block":
            keyword_ = line_s[1].lower()
            #print keyword_
            
            if keyword_ == "latticeparameters":
                abc_cell_block = lines[i+1].split()

            elif keyword_ == "latticevectors":
                vec_block = lines[i+1:i+4]

            elif keyword_ == "atomiccoordinatesandatomicspecies":
                atoms_block = lines[i+1:i+1+atoms_length]
                #print "atoms_block", atoms_block

            elif keyword_ == "chemicalspecieslabel":
                temp = lines[i+1:i+1+n_of_species]
                for spec in temp:
                    species.append(spec.split()[2])
                #print species
        i +=1

    # cell converting
    va = 0; vb = 0; vc = 0
    if (not abc_cell_block) and vec_block:
        a1, a2, a3 = vec_block[0].split()
        a1 = float(a1); a2 = float(a2); a3 = float(a3)
        b1, b2, b3 = vec_block[1].split()
        b1 = float(b1); b2 = float(b2); b3 = float(b3)
        c1, c2, c3 = vec_block[2].split()
        c1 = float(c1); c2 = float(c2); c3 = float(c3)
        va = np.array([a1, a2, a3])
        vb = np.array([b1, b2, b3])
        vc = np.array([c1, c2, c3])
        if cell_scale == 'Ang':
            va = lattice_constant * va
            vb = lattice_constant * vb
            vc = lattice_constant * vc
        elif cell_scale == 'Bohr':
            va = lattice_constant * bohr2ang * va
            vb = lattice_constant * bohr2ang * vb
            vc = lattice_constant * bohr2ang * vc
        else:
            #print "Can`t find cell scale"
            pass

        #a, b, c, alpha, beta, gamma = convert_xyz2abc(va, vb, vc)
        cell = np.array([va,vb,vc])

    elif abc_cell_block and (not vec_block):
        a, b, c, alpha, beta, gamma = abc_cell_block.split()
        a = float(a); b = float(b); c = float(c)
        alpha = float(alpha); beta = float(beta); gamma = float(gamma)
        cell = [a, b, c, alpha, beta, gamma]

    # atoms
    for atm in atoms_block:

        if len(atm.split()) == 4:
            x, y, z, spec = atm.split()
        elif len(atm.split()) == 5:
            x, y, z, spec, serial = atm.split()

        x = float(x); y = float(y); z = float(z); spec = int(spec)

        if _is_ang_scale:
            pass
        elif _is_bohr_scale:
            x = bohr2ang * x; y = bohr2ang * y; z = bohr2ang * z

        elif _is_scaled_ang_scale:
            #if vec_cell:
            x = lattice_constant*x
            y = lattice_constant*y
            z = lattice_constant*z
            #elif not vec_cell:
            #    print "Can`t guess cell scale and type\n"
#        elif _is_fraction_scale:
        
        group = re.findall('\d', species[spec-1])
        sym = ''.join(re.findall('[a-zA-Z]', species[spec-1]))
        if group:
            group = int(''.join(group))
        else:
            group = None
        
        atom = Atom(sym, (x, y, z), groupid=group)
        atoms.append(atom)

    if cell.shape == (3,3):
        #XYZ.write_xyz(file_name.replace('fdf','xyz'), atoms, cell)
        return AtomsSystem(atoms, cell=cell)
    else:
        #XYZ.write_xyz(file_name.replace('fdf','xyz'), atoms)
        return AtomsSystem(atoms, cell=None)

def convert_abc2xyz(a,b,c,alpha,beta,gamma):
    """Convert to cartesian lattice vectors.
       Taken from the routine 
       /biodesign/v330/common/code/source/xtlgraf_batch/celori.f"""
    s1 = sin(alpha*degrad)
    s2 = sin(beta*degrad)
    s3 = sin(gamma*degrad)
    c1 = cos(alpha*degrad)
    c2 = cos(beta*degrad)
    c3 = cos(gamma*degrad)
    c3bar = (-c1*c2 + c3)/(s1*s2)
    sqrtarg = 1.0 - c3bar*c3bar
    if (sqrtarg <= 0):
        print ("Negative argument to SQRT")
        print (sqrtarg, alpha, beta, gamma)
        sqrtarg = 0.0
    s3bar = sqrt(sqrtarg)
    # The general transformation from scaled to XYZ coordinates is now:
    # x = or1 * xabc
    # y = or2 * xabc + or3 * yabc
    # z = or4 * xabc + or5 * yabc + or6 * zabc
    or1 = a*s2*s3bar
    or2 = a*s2*c3bar
    or3 = b*s1
    or4 = a*c2
    or5 = b*c1
    or6 = c
    # Compute Cartesian vectors for a, b, and c
    va = or1, or2, or4
    vb = 0, or3, or5
    vc = 0, 0, or6
    return va,vb,vc

def convert_xyz2abc(va, vb, vc):
    va = Vector(va); vb = Vector(vb); vc = Vector(vc)
    a = va.length(); b = vb.length(); c = vc.length()
    alpha = vb.angle(vc) / pi * 180.
    beta = vc.angle(va) / pi * 180.
    gamma = va.angle(vb) / pi * 180.
    return [a, b, c, alpha, beta, gamma]

def cleansymb(s):
    """
    This function strips off everything after and including the first
    non-letter in an element name.
    Element name is capitalized.
    """
    return (re.split('[^a-zA-Z]',s)[0]).capitalize()

def get_unique_symbs(atoms):
    """
    Get unique symbols from "atoms".
    *IN
    - atoms: list of "atom"s in the format "(Symbol,x,y,z)"
    *OUT
    - unique_symbs
    """
    unique_symbs = []
    for atom in atoms:
        if atom.get_symbol() not in unique_symbs:
            unique_symbs.append(atom.get_symbol())
    return unique_symbs

def read_xyz(file_name, keeptype=False, clean=True, initial=False):
    #090217 <read_ixyz> included ("initial" option) 
    """
    Read the (default=last) structure from an XYZ file.

    Usage:
    >>> atoms = read_xyz(file_name)

    Parameters:
    file_name : the name of a xyz-format file
    Cell inforamtion is available by adding 6 values on comment line (2nd line 
    in a xyz file)
    """
    file_name = file_name.strip()
    try:
        xyz_file = open(file_name,'r')
    except IOError as exception:
        print (exception)
        #...Alternative...
        #type, message, traceback = sys.exc_info()
        #print 'exception type:',type
        #print 'exception message',message
        return
        #sys.exit()
    # First, read the initial structure...
    # line 1: number of atoms
    line = xyz_file.readline()
    words = line.split()
    natm = int(words[0])

    # line 2: title (cell info as custom def.)
    line = xyz_file.readline()
    words = line.split()
    if words:
        if words[0].upper() == 'CELL':
            #print 'This is a CRYSTAL system. (Cell information available)'
            cell = list(map(float,words[1:7]))
        else:
            #print 'This is a MOLECULAR system. (NO cell information)'
            cell = None
    else:
        #print 'This is a MOLECULAR system. (NO cell information)'
        cell = None

    # line 3-: atom, position
    atoms = []
    for i in range(natm):
        line = xyz_file.readline()
        words = line.split()
        symb = words[0]
        if(clean == True):
            symb2 = cleansymb(symb)
        x,y,z = float(words[1]),float(words[2]),float(words[3])
        at = Atom(symb2,[x,y,z], serial=i+1)
        #print symb2
        if(keeptype == True):#111122
            type = symb#111122
            #print keeptype
            #print type
            at.set_fftype(type)#111122
        atoms.append(at)
            
    # If only the initial configuration is needed...
    if initial is True:
        #print 'Read the (initial) structure from "%s" ...' % file_name
        xyz_file.close()
        return AtomsSystem(atoms)

    # Otherwise, repeat until the last structure...
    while True:
        # line 1: number of atoms
        line = xyz_file.readline()
        if not line: break
        words = line.split()
        if not words: break
        natm = int(words[0])
        
        # line 2: title (cell info as custom def.)
        line = xyz_file.readline()
        words = line.split()
        if words[0].upper() == 'CELL':
            cell = list(map(float,words[1:7]))
        else:
            cell = None
            
        # line 3-: atom, position
        atoms = []
        for i in range(natm):
            line = xyz_file.readline()
            words = line.split()
            symb = words[0]
            if clean == True:
                symb = cleansymb(symb)
            x,y,z = float(words[1]),float(words[2]),float(words[3])
            #atoms.append(Atom(symb, [x,y,z]))
            atoms.append(Atom(symb, [x,y,z], serial=i+1))
    #print 'Read the (last) structure from "%s" ...' % file_name            
    xyz_file.close()
    if cell:
        return AtomsSystem(atoms,
                           cell=convert_abc2xyz(cell[0], cell[1], cell[2], 
                                                cell[3], cell[4], cell[5]))
    else:
        return AtomsSystem(atoms)


def write_xyz(file_name, atoms, comm=None, append=False):
    """
    Writes an XYZ file.
    Atom format: (Symbol,x,y,z)
    Cell format: (a,b,c,alpha,beta,gamma)
    """
    file_name = file_name.strip()
    #print 'Writing ..."%s"' % file_name
    if append:
        xyz_file = open(file_name,'a')            
    else:
        xyz_file = open(file_name,'w')
    
    # line 1: number of atoms
    xyz_file.write("%d \n" % len(atoms))

    # line 2: title
    cell = atoms.get_cell()
    if cell is not None:
        a,b,c,alpha,beta,gamma = convert_xyz2abc(cell[0], cell[1], cell[2])
        if comm is not None:
            xyz_file.write("%4s %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %s\n" %\
                           ("CELL",a,b,c,alpha,beta,gamma,comm))
        else:
            xyz_file.write("%4s %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n" %\
                           ("CELL",a,b,c,alpha,beta,gamma))
    else:
        if comm is not None:
            xyz_file.write("%s\n" % comm)
        else:
            xyz_file.write("Insert title information here.\n")
        
    # line 3-: atoms information
    i = 0
    while i < len(atoms._atoms):
        symb = atoms[i].get_symbol()
        x,y,z = atoms[i].get_position()
        #symb,x,y,z = atm
        #xyz_file.write(" %2s %13.6f %13.6f %13.6f\n" %\
        #               (symb,x,y,z))
        xyz_file.write(" %2s %13.6f %13.6f %13.6f\n" %\
                       (symb,x,y,z))
        i += 1
    xyz_file.close()
    return

def parse_lat_line(line):
    words = line.split()
    a,b,c = float(words[1]),float(words[2]),float(words[3])
    alpha,beta,gamma = float(words[4]),float(words[5]),float(words[6])
    return convert_abc2xyz(a,b,c,alpha,beta,gamma)

def parse_atm_line(line):
    words = line.split()
    symb = cleansymb(words[2])
    x,y,z = float(words[6]),float(words[7]),float(words[8])
    fftype = words[9]
    charge = float(words[12])
    return symb,x,y,z, fftype, charge

def parse_con_line(line):
    words = line.split()
    con = []; i = 1
    while i < len(words):
        con.append(int(words[i]))
        i += 1
    return con

def parse_ord_line(line):
    words = line.split()
    ord = []; i = 1
    while i < len(words):
        ord.append(words[i])
        i += 1
    return ord

def read_bgf(file_name):

    try:
        file = open(file_name,'r')
    except IOError as exception:
        print (exception)
        #...Alternative...
        #type, message, traceback = sys.exc_info()
        #print 'exception type:',type
        #print 'exception message',message
        return
    else:
        print ('Reading "%s" ...' % file_name)
   
    latpat = re.compile('CRYSTX')
    atmpat = re.compile('HETATM')
    conpat = re.compile('CONECT')
    ordpat = re.compile('ORDER')

    symbols = []
    atoms = []
    fftypes = []
    charges = []
    cell = []
    connectivities = []
    orders = []
    
    while True:
        line = file.readline()
        if not line: break
        
        if latpat.search(line):
            for v in parse_lat_line(line):
                cell.append(v)

        elif atmpat.search(line):
            symb,x,y,z,fftype,charge = parse_atm_line(line)
            symbols.append(symb)
            atoms.append([x,y,z])
            fftypes.append(fftype)
            charges.append(charge)

        elif conpat.search(line):
            if 'FORMAT' not in line:
                connectivity = parse_con_line(line)
                connectivities.append(connectivity)

        elif ordpat.search(line):
            order = parse_ord_line(line)
            orders.append(order)
    file.close()

    atoms2 = []; i=0
    while i < len(atoms):
        symbol = symbols[i]; position = atoms[i]
        fftype = None; charge = None; connectivity = None
        if fftypes: fftype = fftypes[i]
        if charges: charge = charges[i]
        if connectivities: connectivity = connectivities[i]
        atom = Atom(symbol, position, fftype=fftype, serial=i+1,
                    charge=charge, connectivity=connectivity)
        atoms2.append(atom); i+=1
    if not cell: return AtomsSystem(atoms2)
    else: return AtomsSystem(atoms2, cell=cell)

def write_bgf(file_name, atoms):
    """
    Read the structure from an bgf file.

    Usage:
    >>> atoms = read_bgf(file_name)

    Parameters:
    file_name : the name of a bgf-format file
    """
    #YHK030620
    #print 'Writing "%s" ...' % file_name
    bgf_file = open(file_name,'w')
    cell = atoms.get_cell()
    
    if cell.all():
        cell_ = convert_xyz2abc(cell)
        bgf_file.write("BGFGRF 200\n")
        bgf_file.write("PERIOD 111\n")
        bgf_file.write("AXES ZYX\n")
        bgf_file.write("SGNAME P 1\n")
        bgf_file.write("CRYSTX  %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n"
                       % tuple(cell_))
    else:
        bgf_file.write("BIOGRF 200\n")
        
    bgf_file.write("FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,")
    bgf_file.write("1x,a5,3f10.5,1x,a5, i3,i2,1x,f8.5)\n")

    # atoms
    i = 0
    for atom in atoms._atoms:
        i = i+1
        sym = atom.get_symbol()
        x,y,z = atom.get_position()
        ff_type = sym
        if atom.get_fftype(): ff_type = atom.get_fftype()
        q1 = q2 = 0 # ?????
        chg = 0.
        if atom.get_charge(): chg = atom.get_charge()
        bgf_file.write("%6s %5d %-5s %3s %1s %4s %10.5f%10.5f%10.5f " %
                   ("HETATM",i,sym+str(i),"RES","A","444",x,y,z))#110615
        bgf_file.write("%-5s%3d%2d %8.5f\n" % (ff_type,q1,q2,chg))

    # bondings
    bgf_file.write("FORMAT CONECT (a6,12i6)\n")
    i = 1
    for atom in atoms._atoms:
        conect = atom.get_connectivity()
        if conect:
            line1 = "CONECT"; line2 = "ORDER "
            for a_conect in conect:
                line1 = line1 + "%6s" % str(a_conect)
                line2 = line2 + "%6s" % "1" # bond order (X)
            line1 = line1 + "\n"; line2 = line2 + '\n'
            bgf_file.write(line1)
            if len(conect) > 1: bgf_file.write(line2)
        else: bgf_file.write("CONECT %6s\n" % str(i))
        i += 1
        
    bgf_file.write("END\n")                   
    bgf_file.close()

def read_xsf(file_name):
    """
    Read the structure from an xsf file.

    Usage:
    >>> atoms = read_xsf(file_name)

    Parameters:
    file_name : the name of a xsf-format file
    """
    f = open(file_name)
    lines = f.readlines()
    atoms = []; cell = []
    i_line = 0
    for line in lines:
        keyword = line.split()[0]
        #print keyword
        # crystal
        if keyword == 'CRYSTAL':
            if lines[i_line+1].split()[0] == 'PRIMVEC':
                vec1 = [float(lines[i_line+2].split()[0]),
                        float(lines[i_line+2].split()[1]),
                        float(lines[i_line+2].split()[2])]
                vec2 = [float(lines[i_line+3].split()[0]),
                        float(lines[i_line+3].split()[1]),
                        float(lines[i_line+3].split()[2])]
                vec3 = [float(lines[i_line+4].split()[0]),
                        float(lines[i_line+4].split()[1]),
                        float(lines[i_line+4].split()[2])]
                cell.append(vec1);cell.append(vec2);cell.append(vec3)
                #print cell
        elif keyword == 'PRIMCOORD':
            len_at = int(lines[i_line+1].split()[0])
            atom_block = lines[i_line+2:i_line+2+len_at+1]
            for atom_line in atom_block:
                #print atom_line
                n,x,y,z = atom_line.split()[:4]
                atoms.append(Atom(atomic_symbol[int(n)],
                                  [float(x),float(y),float(z)]))
        i_line += 1
    if cell: atoms = AtomsSystem(atoms, cell=cell)
    else:    atoms = AtomsSystem(atoms)
    return atoms
                

def write_xsf(file_name, atoms):
    """
    For the given AtomsSystem instance, write an .xsf file.

    Usage:
    >>> write_xsf(file_name, atoms)
    """

    xsf = open(file_name, 'w')

    # Molecular structure...
    if atoms.get_cell() is None:
        xsf.write(' ATOMS\n')
        i = 0
        for atom in atoms._atoms:
            x,y,z = atom.get_position()
            info = atomic_number(atom.get_symbol()), x, y, z
            xsf.write(' %2d %12.6f %12.6f %12.6f\n' % info)
            i += 1
        xsf.close()
        
    else:
        xsf.write(' CRYSTAL\n')
        # Primitive lattice vectors (in Angstroms)
        xsf.write(' PRIMVEC\n')
        va,vb,vc = atoms.get_cell()[0],atoms.get_cell()[1],atoms.get_cell()[2] 
        #va, vb, vc = convert_abc2xyz(cell[0], cell[1], cell[2], 
        #                             cell[3], cell[4], cell[5])
        xsf.write(' %12.6f %12.6f %12.6f\n' % tuple(va))
        xsf.write(' %12.6f %12.6f %12.6f\n' % tuple(vb))
        xsf.write(' %12.6f %12.6f %12.6f\n' % tuple(vc))
        # Conventional lattice vectors (in Angstroms)        
        xsf.write(' CONVVEC\n')
        xsf.write(' %12.6f %12.6f %12.6f\n' % tuple(va))
        xsf.write(' %12.6f %12.6f %12.6f\n' % tuple(vb))
        xsf.write(' %12.6f %12.6f %12.6f\n' % tuple(vc))
        # Atomic coord. in a primitive unit cell (in Angstroms)
        xsf.write(' PRIMCOORD\n')
        xsf.write(' %-7d1\n' % len(atoms))
        i = 0
        for atom in atoms._atoms:
            x,y,z = atom.get_position()
            info = atomic_number(atom.get_symbol()), x, y, z
            xsf.write('% 3d %12.6f %12.6f %12.6f\n' % info)
            i += 1
        xsf.close()

def read_pdb(file_name):
    f = open(file_name)
    lines = f.readlines()
    atoms = []; S = []

    for line in lines:
        keyword = line.split()[0]
        args    = line.split()[1:]
        if   keyword == 'CRYST1': pass
        elif keyword == 'SCALE1' or keyword == 'SCALE2' or keyword == 'SCALE3':
            S.append([float(args[0]),float(args[1]),float(args[2])])
        elif keyword == 'HETATM':
            n,sym1,sym2,resi_seq,x,y,z,occ,temp_fac,sym3 = args[:11]
            atoms.append(Atom(sym1, [float(x),float(y),float(z)]))

    cell_m = (np.matrix(S).T)**-1
    atoms = AtomsSystem(atoms, cell=np.array(cell_m))
    return atoms

        
def write_pdb(file_name, atoms):
    """
    xyz --> pdb with cell infomation
    """
    # Cell infomation
    #atoms = read_xyz(file_name)
    cell = atoms.get_cell()
    cell_m = np.matrix(cell)
    S = (cell_m**-1).T
    cell_6 = convert_xyz2abc(*cell)
    List_S=[]
    for element in S.flat:
        List_S.append(element)

    #pdb_file = open(file_name.replace('.xyz','.pdb'), 'w')
    pdb_file = open(file_name, 'w')

    # Header
    pdb_file.write('CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1\n' % tuple(cell_6))
    pdb_file.write('SCALE1    %10.5f%10.5f%10.5f%15.5f\n' % (List_S[0],List_S[1],List_S[2],0.0))
    pdb_file.write('SCALE2    %10.5f%10.5f%10.5f%15.5f\n' % (List_S[3],List_S[4],List_S[5],0.0))
    pdb_file.write('SCALE3    %10.5f%10.5f%10.5f%15.5f\n' % (List_S[6],List_S[7],List_S[8],0.0))

    # Atom lines
    i=1 #Atom Serial Number
    occ=1.00 # Occupancy
    resi_seq=0 #Residue sequence number
    temp_fac=0.00 # Temperature Factor

    for atom in atoms:
        symb = atom.get_symbol()
        x,y,z = atom.get_position()
        line = "HETATM%5i %4s %3s  %4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n" % (i,symb,symb,resi_seq,x,y,z,occ,temp_fac,symb)
        pdb_file.write(line)
        i+=1
    pdb_file.write('END\n')
    pdb_file.close()


def read_cif(file_name):

    # all lines
    lines = open(file_name).readlines()

    keys = []

    # loop 1: get n_keys
    for line in lines:
        words = line.split('_')
        if len(words) <= 1: continue
        if words[1].lower() == 'atom':
            if words[2].lower() == 'site':
                keys.append( '_'.join(words[3:]).replace('\n','') )

    print (keys)

    # number of attributes
    n_keys = len(keys)

    # find index
    i = 0; i_symb = 0; i_x = 0; i_y = 0; i_z = 0
    for key in keys:
        if key == 'type_symbol': i_symb = i
        if key == 'fract_x': i_x = i
        if key == 'fract_y': i_y = i
        if key == 'fract_z': i_z = i
        i += 1

    print (n_keys, i_symb, i_x, i_y, i_z)

    a = 0; b = 0; c = 0; alpha = 0; beta = 0; gamma = 0
    atoms = []

    # loop 2: get cell info
    for line in lines:
        words = line.split('_')
        words = words[:-1] + words[-1].split()
        if len(words) <= 1: continue
        if words[1].lower() == 'cell':
            print (words)
            if words[2].lower() == 'length':
                if words[3].lower() == 'a': a = float(words[4])
                if words[3].lower() == 'b': b = float(words[4])
                if words[3].lower() == 'c': c = float(words[4])
            if words[2].lower() == 'angle':
                if words[3].lower() == 'alpha': alpha = float(words[4])
                if words[3].lower() == 'beta' : beta  = float(words[4])
                if words[3].lower() == 'gamma': gamma = float(words[4])

    # abc2xyz results v2 along y v1 in xy plane
    #print a,b,c,alpha,beta,gamma
    cell = convert_abc2xyz(a,b,c,alpha,beta,gamma)
    print (cell)

    # loop 3: get atom info
    for line in lines:
        words = line.split()
        #print words
        if len(words) <= 1: continue
        if len(words) == n_keys:
            symb = words[i_symb]
            u = float(words[i_x]); v = float(words[i_y]); w = float(words[i_z])
            print (symb, u, v, w)
            x,y,z = Vector(cell[0])*u + Vector(cell[1])*v + Vector(cell[2])*w
            atoms.append( Atom(symb, [x,y,z]) )

    # loop 4: get connectivity

    atoms = AtomsSystem(atoms, cell=cell)
    return atoms   


def read_axyz(file_name, clean=True):
    """
    Read the structures from an animated XYZ file and return "Trajectory" 
    instance, a collection of "AtomsSystem" instance.

    Usage:
    >>> animated_atoms = read_axyz(file_name)

    Parameters:
    file_name : the name of a xyz-format file
    """
    file_name = file_name.strip()
    try:
        xyz_file = open(file_name,'r')
    except IOError as exception:
        print (exception); return

    # line 3-: atom, position
    atoms_s = []
    
    # Otherwise, repeat until the last structure...
    while True:
        # line 1: number of atoms
        line = xyz_file.readline()
        if not line: break
        words = line.split()
        if not words: break
        natm = int(words[0])
        
        # line 2: title (cell info as custom def.)
        line = xyz_file.readline()
        words = line.split()

        #if words[0].split() == []: cell=None
        # 
        #elif words[0].upper() == 'CELL':
        #    cell = map(float,words[1:7])
        #else:
        #    cell = None
        cell = None
            
        # line 3-: atom, position
        atoms = []
        for i in range(natm):
            line = xyz_file.readline()
            words = line.split()
            symb = words[0]
            if clean == True:
                symb = cleansymb(symb)
            x,y,z = float(words[1]),float(words[2]),float(words[3])
            #atoms.append(Atom(symb, [x,y,z]))
            atoms.append(Atom(symb, [x,y,z], serial=i+1))
            
        if cell:
            atoms_s.append(AtomsSystem(atoms,
                                       cell=convert_abc2xyz(cell[0], cell[1],
                                                            cell[2], cell[3],
                                                            cell[4], cell[5])))
        else:
            #print AtomsSystem(atoms)
            #print isinstance(AtomsSystem(atoms), AtomsSystem)
            atoms_s.append(AtomsSystem(atoms))

    xyz_file.close()
    return Trajectory(atoms_s)


def read_axsf(file_name):

    f = open(file_name)
    lines = f.readlines()
    nstep = lines[0]
    lines = lines[1:]
    nstep = int(nstep.split()[-1])
    print ("nstep:", nstep)

    ani = []; cell = []
    i_line = 0
    for line in lines:
        keyword = line.split()[0]
        #print keyword

        if keyword == 'CRYSTAL':
            if lines[i_line+1].split()[0] == 'PRIMVEC':
                vec1 = [float(lines[i_line+2].split()[0]),
                        float(lines[i_line+2].split()[1]),
                        float(lines[i_line+2].split()[2])]
                vec2 = [float(lines[i_line+3].split()[0]),
                        float(lines[i_line+3].split()[1]),
                        float(lines[i_line+3].split()[2])]
                vec3 = [float(lines[i_line+4].split()[0]),
                        float(lines[i_line+4].split()[1]),
                        float(lines[i_line+4].split()[2])]
                cell.append(vec1);cell.append(vec2);cell.append(vec3)
                print ("cell\n", cell)

        elif keyword == 'PRIMCOORD':
            i_step = int(line.split()[1]); print ("STEP:", i_step)
            len_at = int(lines[i_line+1].split()[0])
            atom_block = lines[i_line+2:i_line+2+len_at]
            print ("atom_block\n", atom_block)

            atoms = []
            for atom_line in atom_block:
                print (atom_line)
                n,x,y,z = atom_line.split()[:4]
                atoms.append(Atom(atomic_symbol[int(n)],
                                  [float(x),float(y),float(z)]))
            ani.append(atoms)
            temp = AtomsSystem(atoms); print (temp)

        i_line += 1

    ani2 = []
    for at in ani:
        ani2.append(AtomsSystem(at, cell=cell))

    return ani2


def write_axsf(file_name, ani):
    """
    For the given AtomsSystem instance, write an .xsf file.

    Usage:
    >>> write_axsf(file_name, ani)
    """

    # read file
    xsf = open(file_name, 'w')
    xsf.write('ANIMSTEPS %i\n' % len(ani))

    # crystal?
    atoms = ani[0].copy()

    if atoms.get_cell() is None:
        pass

    else: 
        xsf.write('CRYSTAL\n')
        # Primitive lattice vectors (in Angstroms)
        xsf.write('PRIMVEC\n')
        va,vb,vc = atoms.get_cell()[0],atoms.get_cell()[1],atoms.get_cell()[2] 
        xsf.write('%12.6f %12.6f %12.6f\n' % tuple(va))
        xsf.write('%12.6f %12.6f %12.6f\n' % tuple(vb))
        xsf.write('%12.6f %12.6f %12.6f\n' % tuple(vc))

    # write coord.
    i_step = 0
    while i_step < len(ani):
        atoms = ani[i_step].copy()

        if atoms.get_cell() is None:
            xsf.write('ATOMS %i\n' % (i_step+1))

        else:
            xsf.write('PRIMCOORD %i\n' % (i_step+1))
            xsf.write('%s 1\n' % str(len(atoms)))
        i = 0
        for atom in atoms._atoms:
            x,y,z = atom.get_position()
            info = atomic_number(atom.get_symbol()), x, y, z
            xsf.write('%2d %12.6f %12.6f %12.6f\n' % info)
            i += 1
        i_step += 1

    xsf.close()

def read_xdat(file_name):
    #path=os.getcwd()+'/XDATCAR'
    #f = open(path)
    f = open(file_name)
    lines = f.readlines()

    # system info.
    line_title = lines[0]
    line_cell1 = lines[2].split()
    line_cell2 = lines[3].split()
    line_cell3 = lines[4].split()
    line_symb  = lines[5].split()
    line_numb  = lines[6].split()

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
        cell1.append(float(v1))
    for v2 in line_cell2:
        cell2.append(float(v2))
    for v3 in line_cell3:
        cell3.append(float(v3))
    cell = [cell1, cell2, cell3]

    # atoms info.
    line_atoms = lines[7:]; snapshots = []
    i = 0
    for line in line_atoms:
        atoms = []
        i += 1
        if line.split()[0].lower() == 'direct':
            n_step = int(line.split()[-1])
            line_coord = line_atoms[i:i+n_system]
            j = 0
            for coord in line_coord:
                xf,yf,zf = coord.split()
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
            #name = 'conf_%6.6i.xyz' % n_step
            #io.write_xyz(name, atoms_obj)
            ### snapshots => "list" of AtomsSystem instances
            snapshots.append(atoms_obj)
    return snapshots
    
