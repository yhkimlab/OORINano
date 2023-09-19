"""
2023.06.26: updating following refactoring: 
"""

from ..atoms import *
from ..ncio import cleansymb, get_unique_symbs
from ..units import ang2bohr; bohr2ang = 1/ang2bohr


def make_seqquest_positions(atoms):
    # positions
    atoms1 = []
    atoms1.append("atom types\n")
    unique_symbs = get_unique_symbs(atoms)
    atoms1.append(" %d \n" % len(unique_symbs))
    for symb in unique_symbs:
        atoms1.append("atom file\n")
        atoms1.append(" %s.atm\n" % symb)
    atoms1.append("number of atoms in unit cell\n")
    atoms1.append(" %d\n" % len(atoms))
    atoms1.append("atom, type, position vector\n")
    n = 1
    for atom in atoms:
        symb = atom.get_symbol(); x,y,z = atom.get_position()
        atoms1.append(" %4d %3d %15.9f %15.9f %15.9f\n" % \
                     (n,unique_symbs.index(symb)+1,
                      x*ang2bohr,y*ang2bohr,z*ang2bohr))
        n = n + 1
    return atoms1

def write_seqquest_input(atoms, ref_file):
    """
    make 'lcao.in' using a reference file
    IN
    - cellxyz : a list from xyz2quest.convert_abc2xyz()
    - atoms : a list generated by XYZ.read_xyz()
    - ref_file : use its contents except structure part
    OUT
    - converted.in 
    """
    # positions
    atoms1 = make_seqquest_positions(atoms)

    # parameters from a reference file
    ref_lines1 = []; ref_lines2 = []
    ref = open(ref_file)
    ref_lines = ref.readlines()
    for n in range(len(ref_lines)):
        if ref_lines[n] == 'grid dimensions\n': # ref_lines1
            ref_lines1 = ref_lines[:(n+2)]
        elif ref_lines[n] == 'kgrid\n':  # ref_lines2
            ref_lines2 = ref_lines[n:]; break
        elif ref_lines[n] == 'end setup phase data':
            ref_lines2 = ref_lines[n:]; break
        else:
            pass
            
    file = open('converted.in','w')
    for x in ref_lines1, atoms1, ref_lines2:
        for y in x:
            file.write(y)
    file.close()

# DFT parameters
par_op = {'force':1,                   # (true or false)
          'atom_relax':1,              # (true or false)
          'cell_relax':0,              # (true or false)
          'neb':0,                     # (true or false)
          'title':'Insert Title Here',
          'XC':'PBE',                  # LDA, LDA-SP, PBE, PBE-SP, PW91, BLYP
          #'spin_polarization':0.,      # units of excess e- of majoruty spin
          #'Efield':[0.,0.,0.],         # Ext. E field vector [Ry/bohr]
          #'epsilon':0.,                # bulk material static dielectric
          #'strain':[[11,21,31],[],[]], # deformation matrix strain(i,j)
          'kgrid':[1,1,0],             # k point sampling
          }
par_en = {'force':0,                   # (true or false)
          'atom_relax':0,              # (true or false)
          'cell_relax':0,              # (true or false)
          'neb':0,                     # (true or false)
          'title':'Insert Title Here',
          'XC':'PBE',                  # LDA, LDA-SP, PBE, PBE-SP, PW91, BLYP
          #'spin_polarization':0.,      # units of excess e- of majoruty spin
          #'Efield':[0.,0.,0.],         # Ext. E field vector [Ry/bohr]
          #'epsilon':0.,                # bulk material static dielectric
          #'strain':[[11,21,31],[],[]], # deformation matrix strain(i,j)
          'kgrid':[1,1,0],             # k point sampling
          }


def write_seqquest(atoms, params, auto_sort=1, file_name=''):
    """
    make lcao.in according to current model and parameters


    usage:
    
    >>> write_seqquest(atoms, params, auto_sort=1, file_name='')

    * IN
    - atoms : AtomsSystem instance
    - params : parameters for SeqQuest (refer to quest.par_op/quest.par_en)
    - auto_sort : whether atomic coordinations will be sorted inside the
      simulation box or not, due to the specific convention of SeqQuest
      (default=True)
    - file_name (default='lcao.in_by_XXYZ')

    * OUT
    - No return value (file output)
    """
    from math import sqrt
    if not file_name:
        file = open('lcao.in_by_XXYZ','w')
    else:
        if type(file_name) == str:
            file = open(file_name, 'w')
        else:
            raise ValueError ("file_name : string type is expected.")   # error by pdoc, J Park
    # commends
    file.write("do setup\n")
    file.write("do iters\n")
    file.write("do post\n")
    if params['force']: file.write("do force\n") # force calculation?
    else:               file.write("no force\n")
    if params['atom_relax']: file.write("do relax\n") # geometry relaxation?
    else:                    file.write("no relax\n")
    if params['cell_relax']: file.write("do cell\n") # cell relaxation?
    else:                    file.write("no cell\n")
    if params['neb']: file.write("do neb\n") # NEB transition state finding
    else:             file.write("no neb\n")

    # setup phase data
    file.write("setup data\n")
    file.write("title\n")
    file.write(" %s\n" % params['title'])
    file.write("functional\n")
    file.write(" %s\n" % params['XC']) # exchange-correlation functional
    file.write("dimension\n")
    pbc_x, pbc_y, pbc_z = atoms.get_pbc()
    atoms2 = atoms.copy(); atoms2.select_all()
    atoms_org = atoms.copy(); atoms_org.select_all()

    # modified scheme (110629)
    if auto_sort:
        vec1, vec2, vec3 = atoms2.get_cell()
        atoms2.sort('x'); avx = (atoms2[0][0]+atoms2[-1][0])/2
        atoms2.sort('y'); avy = (atoms2[0][1]+atoms2[-1][1])/2
        atoms2.sort('z'); avz = (atoms2[0][2]+atoms2[-1][2])/2
        atoms2.translate(-avx,-avy,-avz)
        atoms_org.translate(-avx,-avy,-avz)

        if atoms.get_pbc().all():
            file.write(" 3\n")
            vec_t = 0.5*(vec1+vec2+vec3)
            atoms2.translate(*vec_t)
            atoms_org.translate(*vec_t)
        elif pbc_x and pbc_y:
            file.write(" 2\n")
            vec_t = 0.5*(vec1+vec2)
            atoms2.translate(*vec_t)
            atoms_org.translate(*vec_t)
        elif pbc_x:
            file.write(" 1\n")
            vec_t = 0.5*vec1
            atoms2.translate(*vec_t)
            atoms_org.translate(*vec_t)
        else:
            file.write(" 0\n") #warning!
    else:
         if atoms.get_pbc().all(): file.write(" 3\n")
         elif pbc_x and pbc_y:     file.write(" 2\n")
         elif pbc_x:               file.write(" 1\n")
         else:                     file.write(" 0\n")
         
    file.write("primitive lattice vectors\n")
    va = Vector([20.0, 0.0, 0.0])
    vb = Vector([0.0, 20.0, 0.0])
    vc = Vector([0.0, 0.0, 20.0]) #warning!
    if atoms.get_cell() is not None:
        va,vb,vc = atoms.get_cell()
        va = Vector(va); vb = Vector(vb); vc = Vector(vc)
    file.write("%15.9f %15.9f %15.9f\n" % tuple(ang2bohr*va))
    file.write("%15.9f %15.9f %15.9f\n" % tuple(ang2bohr*vb))
    file.write("%15.9f %15.9f %15.9f\n" % tuple(ang2bohr*vc))
    file.write("grid dimensions\n")
    h = 0.3 # grid size
    n1 = int( ang2bohr*sqrt(va[0]**2+va[1]**2+va[2]**2)/h ) + 1
    n2 = int( ang2bohr*sqrt(vb[0]**2+vb[1]**2+vb[2]**2)/h ) + 1
    n3 = int( ang2bohr*sqrt(vc[0]**2+vc[1]**2+vc[2]**2)/h ) + 1
    file.write(" %i %i %i\n" % (n1,n2,n3))

    # positions
    atoms1 = make_seqquest_positions(atoms_org)
    for line in atoms1: file.write(line)
    if atoms.get_pbc().any():
        file.write("kgrid\n")
        file.write(" %i %i %i\n" % tuple(params['kgrid'])) # k point sampling
    file.write("end setup phase data\n")

    # run phase data
    file.write("run phase data\n")
    file.write("convergence\n")
    file.write(" 0.000500\n")
    #file.write("history\n")
    #file.write(" 10\n")

    # geometry
    if params['atom_relax']:
        file.write("geometry parameters\n")
        file.write("gmethod\n")
        file.write(" ASD\n")
        file.write("gsteps\n")
        file.write(" 50\n")
        file.write("gconv\n")
        file.write(" 0.0010\n")
        file.write("end geometry parameters\n")

    # cell optimization
    if params['cell_relax']:
        file.write("ucsteps\n")
        file.write(" 50\n")
        file.write("uc_convergence\n")
        file.write("0.0010\n")
        file.write("ucmethod\n")
        file.write(" STEEPEST\n")
        
    file.write("end run phase data\n")


def read_seqquest(file_name='lcao.in', get_params = 0):
    """
    get information about atomic coordinations and DFT parameters from lcao.in


    usage:

    >>> atoms = read_seqquest()
    >>> atoms = read_seqquest('lcao.in_my')
    >>> atoms = read_seqquest('lcao.in_my', get_params=0)
    >>> atoms, par = read_seqquest(file_name='lcao.in_my', get_params=1)

    *IN
    - file_name : input file of SeqQuest (default='lcao.in')
    - get_params : whether you want to get DFT parameters or not, (default=False)

    *OUT
    - atoms : AtomsSystem instance
    - par : DFT parameters for SeqQuest (refer to quest.par_op/quest.par_en)
    """
    
    param = {'force':0,                  # (true or false)
             'atom_relax':0,             # (true or false)
             'cell_relax':0,             # (true or false)
             'neb':0,                    # (true or false)
             'title':'Insert Title Here',
             'XC':'PBE',                 # LDA, LDA-SP, PBE, PBE-SP, PW91, BLYP
             #'spin_polarization':0.,     # units of excess e- of majoruty spin
             #'Efield':[0.,0.,0.],        # Ext. E field vector [Ry/bohr]
             #'epsilon':0.,               # bulk material static dielectric
             #'strain':[[11,21,31],[],[]],# deformation matrix strain(i,j)
             'kgrid':[1,1,0],            # k point sampling
             }

    # parameters
    dim = 0
    v1 = 0; v2 = 0; v3 = 0
    n_species = 0; i_species = 1; species = {}
    n_atoms = 0; atoms = []

    # read lcao.in
    f = open(file_name)
    lines = f.readlines()

    i_line = 0
    for line in lines:
        ## DFT parameters
        if line == 'do force\n':    param['force'] = 1
        elif line == 'no force\n':  param['force'] = 0
        elif line == 'do relax\n':  param['atom_relax'] = 1
        elif line == 'no relax\n':  param['atom_relax'] = 0
        elif line == 'do cell\n':   param['cell_relax'] = 1
        elif line == 'no cell\n':   param['cell_relax'] = 0
        elif line == 'do neb\n':    param['neb'] = 1
        elif line == 'no neb\n':    param['neb'] = 0
        elif line == 'title\n':     param['title'] = lines[i_line+1]
        elif line == 'functional\n':param['XC'] = lines[i_line+1].split()[0]
        elif line == 'kgrid\n':
            kx,ky,kz = lines[i_line+1].split()
            param['kgrid'] = [int(kx), int(ky), int(kz)]
            
        ## AtomsSystem
        # dimension
        elif line == 'dimension\n':
            dim = int(lines[i_line+1].split()[0])
        # cell info.
        elif line == 'primitive lattice vectors\n':
            v1 = lines[i_line+1].split()
            v1 = bohr2ang*Vector([float(v1[0]), float(v1[1]), float(v1[2])])
            v2 = lines[i_line+2].split()
            v2 = bohr2ang*Vector([float(v2[0]), float(v2[1]), float(v2[2])])
            v3 = lines[i_line+3].split()
            v3 = bohr2ang*Vector([float(v3[0]), float(v3[1]), float(v3[2])])
        # number of species
        elif line == 'atom types\n':
            n_species = int(lines[i_line+1].split()[0])
        # elements
        elif line == 'atom file\n':
            sym = lines[i_line+1].split('.')[0].split()[0]
            species[i_species] = sym
            i_species += 1
        # number of atoms
        elif line == 'number of atoms in unit cell\n':
            n_atoms = int(lines[i_line+1].split()[0])
        # positions
        elif line == 'atom, type, position vector\n':
            atom_lines = lines[i_line+1:i_line+1+n_atoms]
        i_line += 1

    for atline in atom_lines:
        serial,spec,x,y,z = atline.split()
        spec = int(spec)
        x = bohr2ang*float(x); y = bohr2ang*float(y); z = bohr2ang*float(z) 
        atoms.append( Atom(species[spec], (x,y,z)) )

    if not get_params: return AtomsSystem(atoms, cell=[v1,v2,v3], pbc=dim)
    else: return AtomsSystem(atoms, cell=[v1,v2,v3], pbc=dim), param


def select_rngnbs(astr):
    """
    Return the list of selected positive integer numbers.
    *IN 
    - astr: A str. specifying atom #s to be selected. Can be in
            the range form expressed via '-' (e.g. '5-7', '10 - 12')
    *OUT
    - selected: A list of "selected" atom numbers
    (e.g.) '1-3 20 10 - 12 5' => [1, 2, 3, 5, 10, 11, 12, 20]
    """
    # Check whether illegal non-digital characters are present.
    nondgts = re.findall(r'\D+',astr)
    for i in nondgts:
        j = i.strip()
        if j != '' and j != '-':
            raise ValueError('Illegal non-digital character(s)=%s' % j)     # error by pdoc: J Park
    # First, atom #s given in range formats (e.g. '5-7', '10 - 15')
    rngs = re.findall(r'(\d+)\s*-\s*(\d+)', astr)
    selected = []
    for (lower,upper) in rngs:
        num = int(lower)
        while num <= int(upper):
            selected.append(num)
            num += 1
    # Second, atom #s given individually
    all = re.findall(r'\d+',astr)
    for i in all:
        num = int(i)
        if num not in selected: selected.append(num)
    # Finally, sort them
    selected.sort()
    return selected

