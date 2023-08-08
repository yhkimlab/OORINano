''' atomic properties in new object 
atom_prop: atom symbol(key) = [ atomic number, magnetic moment]
    MAGMOM in INCAR: magnetic moment*1.5 is applied
'''

atom_prop = {
    'H':    [1,  1],
    'He':   [2,  0],
    'Li':   [3,  1],
    'Be':   [4,  2],
    'B':    [5,  3],
    'C':    [6,  1],        # valance e = 4, normally 0, if radical 1
    'N':    [7,  1],        # valance e = 3, normally 0, if radical 1
    'O':    [8,  2],
    'F':    [9,  1],
    'Ne':   [10, 0],
    'Na':   [11, 1],
    'Mg':   [12, 2],
    'Al':   [13, 3],
    'Si':   [14, 4],
    'P':    [15, 3],
    'S':    [16, 2],
    'Cl':   [17, 1],
    'Ar':   [18, 0],
    'K':    [19, 1],
    'Ca':   [20, 2],
    'Sc':   [21, 3],
    'Ti':   [22, 4],
    'V':    [23, 5],
    'Cr':   [24, 6],
    'Mn':   [25, 5],
    'Fe':   [26, 4],
    'Co':   [27, 3],
    'Ni':   [28, 2],
    'Cu':   [29, 1],
    'Zn':   [30, 0],
    'Pt':   [78, 2],
    }
