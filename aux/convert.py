'''
    ln -s env_yours.py env.py
        link softly your envfile to env.py
'''

from __future__ import print_function
import re, sys
from math import sin, cos, sqrt, pi
from ..units import ang2bohr, degrad
from ..atoms import atomic_weight, atomic_symbol, atomic_number
from ..atoms.atom_vector import Vector

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

def convert_time2human(sec):
    day = sec // (24 * 3600)
    sec %= 24 * 3600
    hour = sec // 3600
    sec %= 3600
    mini = sec // 60
    sec %= 60
    st = f"{day}d - {hour:d}:{mini:02d}:{sec:02d}"
    return st

def check_file(fname, st):
    with open(fname, 'r') as f:
        for line in f.readlines():
            if st in line:
                return True
    return False
