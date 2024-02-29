'''
auxiliary functions for catalysis
'''

def fix_atoms(atoms, fix):
    if fix == 'b1L':
        atoms.select_atoms("gid", 0)
    else:
        pass
    return None