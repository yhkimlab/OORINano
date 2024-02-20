#
# Object-Oriented XXYZ: atomic data taken from ASE
#
# atomic data : atomic number, symbol, and mass
#               crystal information in reference state
#
# Last revision : 
#
# Revision history
#


atomic_symbol = {
#Period 1   
 1:'H',                                                         2:'He',
#Period 2
 3:'Li',  4:'Be',  5:'B',   6:'C',   7:'N',   8:'O',   9:'F',  10:'Ne',
#Period 3
11:'Na', 12:'Mg', 13:'Al', 14:'Si', 15:'P',  16:'S',  17:'Cl', 18:'Ar',
#Period 4
19:'K',  20:'Ca',
21:'Sc',22:'Ti',23:'V',24:'Cr',25:'Mn',26:'Fe',27:'Co',28:'Ni',29:'Cu',30:'Zn',
                  31:'Ga', 32:'Ge', 33:'As', 34:'Se', 35:'Br', 36:'Kr',
#Period 5
37:'Rb', 38:'Sr',
39:'Y',40:'Zr',41:'Nb',42:'Mo',43:'Tc',44:'Ru',45:'Rh',46:'Pd',47:'Ag',48:'Cd',
                  49:'In', 50:'Sn', 51:'Sb', 52:'Te', 53:'I',  54:'Xe',
#Period 6
55:'Cs', 56:'Ba',
71:'Lu',72:'Hf',73:'Ta',74:'W',75:'Re',76:'Os',77:'Ir',78:'Pt',79:'Au',80:'Hg',
                  81:'Ti', 82:'Pb', 83:'Bi', 84:'Po', 85:'At', 86:'Rn',
#Period 7
87:'Fr', 88:'Ra',
104:'Rf',105:'Db',106:'Sg',107:'Bh',108:'Hs',109:'Mt',110:'Ds',111:'Rg',
#Lanthanoids
57:'La',58:'Ce',59:'Pr',60:'Nd',61:'Pm',62:'Sm',63:'Eu',64:'Gd',65:'Tb',
66:'Dy',67:'Ho',68:'Er',69:'Tm',70:'Yb',
#Actinoids
89:'Ac',90:'Th',91:'Pa',92:'U',93:'Np',94:'Pu',95:'Am',96:'Cm',97:'Bk',98:'Cf',
99:'Es',100:'Fm',101:'Md',102:'No',103:'Lr',

# NMJ MODIFICATION
104:'Had', 105:'Oad'
}

hatomic_number = {
#Period 1   
 'H':1,                                                         'He':2,
#Period 2
 'Li':3,  'Be':4,  'B':5,   'C':6,   'N':7,   'O':8,   'F':9,  'Ne':10,
#Period 3
'Na':11, 'Mg':12, 'Al':13, 'Si':14, 'P':15,  'S':16,  'Cl':17, 'Ar':18,
#Period 4
'K':19,  'Ca':20,
'Sc':21,'Ti':22,'V':23,'Cr':24,'Mn':25,'Fe':26,'Co':27,'Ni':28,'Cu':29,'Zn':30,
                  'Ga':31, 'Ge':32, 'As':33, 'Se':34, 'Br':35, 'Kr':36,
#Period 5
'Rb':37, 'Sr':38,
'Y':39,'Zr':40,'Nb':41,'Mo':42,'Tc':43,'Ru':44,'Rh':45,'Pd':46,'Ag':47,'Cd':48,
                  'In':49, 'Sn':50, 'Sb':51, 'Te':52, 'I':53,  'Xe':54,
#Period 6
'Cs':55, 'Ba':56,
'Lu':71,'Hf':72,'Ta':73,'W':74,'Re':75,'Os':76,'Ir':77,'Pt':78,'Au':79,'Hg':80,
                  'Ti':81, 'Pb':82, 'Bi':83, 'Po':84, 'At':85, 'Rn':86,
#Period 7
'Fr':87, 'Ra':88,
'Rf':104,'Db':105,'Sg':106,'Bh':107,'Hs':108,'Mt':109,'Ds':110,'Rg':111,
#Lanthanoids
'La':57,'Ce':58,'Pr':59,'Nd':60,'Pm':61,'Sm':62,'Eu':63,'Gd':64,'Tb':65,
'Dy':66,'Ho':67,'Er':68,'Tm':69,'Yb':70,
#Actinoids
'Ac':89,'Th':90,'Pa':91,'U':92,'Np':93,'Pu':94,'Am':95,'Cm':96,'Bk':97,'Cf':98,
'Es':99,'Fm':100,'Md':101,'No':102,'Lr':103,

# NMJ MODIFICATION
'Had':104, 'Oad':105
}

def atomic_number(symb):
    if symb == 'X': return 0
    symb = symb.capitalize()
    #for k,v in atomic_symbol.iteritems():
    for k,v in atomic_symbol.items():
        if v==symb: return k
    raise ValueError('"%s" is not defined in "atomic_symbol"' % symb)

# From ASE ver. 3.0
import numpy as N
atomic_weight = N.array([
   0.00000, # X
   1.00794, # H
   4.00260, # He
   6.94100, # Li
   9.01218, # Be
  10.81100, # B
  12.01100, # C
  14.00670, # N
  15.99940, # O
  18.99840, # F
  20.17970, # Ne
  22.98977, # Na
  24.30500, # Mg
  26.98154, # Al
  28.08550, # Si
  30.97376, # P
  32.06600, # S
  35.45270, # Cl
  39.94800, # Ar
  39.09830, # K
  40.07800, # Ca
  44.95590, # Sc
  47.88000, # Ti
  50.94150, # V
  51.99600, # Cr
  54.93800, # Mn
  55.84700, # Fe
  58.93320, # Co
  58.69340, # Ni
  63.54600, # Cu
  65.39000, # Zn
  69.72300, # Ga
  72.61000, # Ge
  74.92160, # As
  78.96000, # Se
  79.90400, # Br
  83.80000, # Kr
  85.46780, # Rb
  87.62000, # Sr
  88.90590, # Y
  91.22400, # Zr
  92.90640, # Nb
  95.94000, # Mo
    N.nan, # Tc
 101.07000, # Ru
 102.90550, # Rh
 106.42000, # Pd
 107.86800, # Ag
 112.41000, # Cd
 114.82000, # In
 118.71000, # Sn
 121.75700, # Sb
 127.60000, # Te
 126.90450, # I
 131.29000, # Xe
 132.90540, # Cs
 137.33000, # Ba
 138.90550, # La
 140.12000, # Ce
 140.90770, # Pr
 144.24000, # Nd
    N.nan, # Pm
 150.36000, # Sm
 151.96500, # Eu
 157.25000, # Gd
 158.92530, # Tb
 162.50000, # Dy
 164.93030, # Ho
 167.26000, # Er
 168.93420, # Tm
 173.04000, # Yb
 174.96700, # Lu
 178.49000, # Hf
 180.94790, # Ta
 183.85000, # W
 186.20700, # Re
 190.20000, # Os
 192.22000, # Ir
 195.08000, # Pt
 196.96650, # Au
 200.59000, # Hg
 204.38300, # Tl
 207.20000, # Pb
 208.98040, # Bi
    N.nan, # Po
    N.nan, # At
    N.nan, # Rn
    N.nan, # Fr
 226.02540, # Ra
    N.nan, # Ac
 232.03810, # Th
 231.03590, # Pa
 238.02900, # U
 237.04820, # Np
    N.nan, # Pu
    N.nan, # Am
    N.nan, # Cm
    N.nan, # Bk
    N.nan, # Cf
    N.nan, # Es
    N.nan, # Fm
    N.nan, # Md
    N.nan, # No
    N.nan, # Lr
    1.00794, # NMJ H for HER
    15.99940])# NMJ O for HER

covalent_radii = N.array([
 0.20, # X
 0.32, # H
 0.93, # He
 1.23, # Li
 0.90, # Be
 0.82, # B
 0.77, # C
 0.75, # N
 0.73, # O
 0.72, # F
 0.71, # Ne
 1.54, # Na
 1.36, # Mg
 1.18, # Al
 1.11, # Si
 1.06, # P
 1.02, # S
 0.99, # Cl
 0.98, # Ar
 2.03, # K
 1.74, # Ca
 1.44, # Sc
 1.32, # Ti
 1.22, # V
 1.18, # Cr
 1.17, # Mn
 1.17, # Fe
 1.16, # Co
 1.15, # Ni
 1.17, # Cu
 1.25, # Zn
 1.26, # Ga
 1.22, # Ge
 1.20, # As
 1.16, # Se
 1.14, # Br
 1.89, # Kr
 2.16, # Rb
 1.91, # Sr
 1.62, # Y
 1.45, # Zr
 1.34, # Nb
 1.30, # Mo
 1.27, # Tc
 1.25, # Ru
 1.25, # Rh
 1.28, # Pd
 1.34, # Ag
 1.41, # Cd
 1.44, # In
 1.41, # Sn
 1.40, # Sb
 1.36, # Te
 1.33, # I
 1.31, # Xe
 2.35, # Cs
 1.98, # Ba
 1.25, # La
 1.65, # Ce
 1.65, # Pr
 1.64, # Nd
 1.63, # Pm
 1.62, # Sm
 1.85, # Eu
 1.61, # Gd
 1.59, # Tb
 1.59, # Dy
 1.58, # Ho
 1.57, # Er
 1.56, # Tm
 1.70, # Yb
 1.56, # Lu
 1.44, # Hf
 1.34, # Ta
 1.30, # W
 1.28, # Re
 1.26, # Os
 1.27, # Ir
 1.30, # Pt
 1.34, # Au
 1.49, # Hg
 1.48, # Tl
 1.47, # Pb
 1.46, # Bi
 1.53, # Po
 1.47, # At
 N.nan, # Rn
 N.nan, # Fr
 N.nan, # Ra
 N.nan, # Ac
 1.65, # Th
 N.nan, # Pa
 1.42, # U
 N.nan, # Np
 N.nan, # Pu
 N.nan, # Am
 N.nan, # Cm
 N.nan, # Bk
 N.nan, # Cf
 N.nan, # Es
 N.nan, # Fm
 N.nan, # Md
 N.nan, # No
 N.nan, # Lr
 0.32,  # NMJ H for HER
 0.73]) # NMJ H for ORR

# This data is from Ashcroft and Mermin.
reference_state = {\
    'X' :{'magmom': 0,                                                                    },                                                                       
    'H' :{'magmom': 1,  'symmetry': 'diatom', 'd': 0.74,                                  },                                          
    'He':{'magmom': 0,  'symmetry': 'atom'                                                },                                                       
    'Li':{'magmom': 1,  'symmetry': 'BCC', 'a': 3.49                                      },                                             
    'Be':{'magmom': 2,  'symmetry': 'hcp', 'c/a': 1.567, 'a': 2.29                        },           
    'B' :{'magmom': 3,  'symmetry': 'Tetragonal', 'c/a': 0.576, 'a': 8.73                 },    
    'C' :{'magmom': 4,  'symmetry': 'Diamond', 'a': 3.57                                  }, # valance e = 4, normally 0, if radical 1
    'N' :{'magmom': 3,  'symmetry': 'diatom', 'd': 1.10                                   }, # valance e = 3, normally 0, if radical 1
    'O' :{'magmom': 2,  'symmetry': 'diatom', 'd': 1.21                                   }, 
    'F' :{'magmom': 1,  'symmetry': 'diatom', 'd': 1.42                                   }, 
    'Ne':{'magmom': 0,  'symmetry': 'fcc', 'a': 4.43                                      },    
    'Na':{'magmom': 1,  'symmetry': 'BCC', 'a': 4.23                                      },    
    'Mg':{'magmom': 2,  'symmetry': 'hcp', 'c/a': 1.624, 'a': 3.21                        },   
    'Al':{'magmom': 3,  'symmetry': 'fcc', 'a': 4.05                                      },                 
    'Si':{'magmom': 4,  'symmetry': 'Diamond', 'a': 5.43                                  },             
    'P ':{'magmom': 3,  'symmetry': 'Cubic', 'a': 7.17                                    },               
    'S' :{'magmom': 2,  'symmetry': 'Orthorhombic', 'c/a': 2.339, 'a': 10.47,'b/a': 1.229 },
    'Cl':{'magmom': 1,  'symmetry': 'Orthorhombic', 'c/a': 1.324, 'a': 6.24, 'b/a': 0.718 },
    'Ar':{'magmom': 0,  'symmetry': 'fcc', 'a': 5.26                                      },
    'K' :{'magmom': 1,  'symmetry': 'BCC', 'a': 5.23                                      },
    'Ca':{'magmom': 2,  'symmetry': 'fcc', 'a': 5.58                                      },
    'Sc':{'magmom': 3,  'symmetry': 'hcp', 'c/a': 1.594, 'a': 3.31                        },
    'Ti':{'magmom': 4,  'symmetry': 'hcp', 'c/a': 1.588, 'a': 2.95                        },
    'V' :{'magmom': 5,  'symmetry': 'BCC', 'a': 3.02                                      },
    'Cr':{'magmom': 6,  'symmetry': 'BCC', 'a': 2.88                                      },
    'Mn':{'magmom': 5,  'symmetry': 'Cubic', 'a': 8.89                                    },
    'Fe':{'magmom': 4,  'symmetry': 'BCC', 'a': 2.87                                      },
    'Co':{'magmom': 3,  'symmetry': 'hcp', 'c/a': 1.622, 'a': 2.51                        },
    'Ni':{'magmom': 2,  'symmetry': 'fcc', 'a': 3.52                                      },
    'Cu':{'magmom': 1,  'symmetry': 'fcc', 'a': 3.61                                      },
    'Zn':{'magmom': 0,  'symmetry': 'hcp', 'c/a': 1.856, 'a': 2.66                        },
    'Ga':{'magmom': 3,  'symmetry': 'Orthorhombic', 'c/a': 1.695, 'a': 4.51, 'b/a': 1.001 },
    'Ge':{'magmom': 4,  'symmetry': 'Diamond', 'a': 5.66                                  },
    'As':{'magmom': 3,  'symmetry': 'Rhombohedral', 'a': 4.13, 'alpha': 54.10             },
    'Se':{'magmom': 2,  'symmetry': 'hcp', 'c/a': 1.136, 'a': 4.36                        },
    'Br':{'magmom': 1,  'symmetry': 'Orthorhombic', 'c/a': 1.307, 'a': 6.67, 'b/a': 0.672 },
    'Kr':{'magmom': 0,  'symmetry': 'fcc', 'a': 5.72                                      },
    'Rb':{'magmom': 1,  'symmetry': 'BCC', 'a': 5.59                                      },
    'Sr':{'magmom': 2,  'symmetry': 'fcc', 'a': 6.08                                      },
    'Y' :{'magmom': 3,  'symmetry': 'hcp', 'c/a': 1.571, 'a': 3.65                        },
    'Zr':{'magmom': 4,  'symmetry': 'hcp', 'c/a': 1.593, 'a': 3.23                        },
    'Nb':{'magmom': 5,  'symmetry': 'BCC', 'a': 3.30                                      },
    'Mo':{'magmom': 6,  'symmetry': 'BCC', 'a': 3.15                                      },
    'Tc':{'magmom': 5,  'symmetry': 'hcp', 'c/a': 1.604, 'a': 2.74                        },
    'Ru':{'magmom': 4,  'symmetry': 'hcp', 'c/a': 1.584, 'a': 2.70                        },
    'Rh':{'magmom': 3,  'symmetry': 'fcc', 'a': 3.80                                      },
    'Pd':{'magmom': 2,  'symmetry': 'fcc', 'a': 3.89                                      },
    'Ag':{'magmom': 1,  'symmetry': 'fcc', 'a': 4.09                                      },
    'Cd':{'magmom': 0,  'symmetry': 'hcp', 'c/a': 1.886, 'a': 2.98                        },
    'In':{'magmom': 3,  'symmetry': 'Tetragonal', 'c/a': 1.076, 'a': 4.59                 },
    'Sn':{'magmom': 4,  'symmetry': 'Tetragonal', 'c/a': 0.546, 'a': 5.82                 },
    'Sb':{'magmom': 3,  'symmetry': 'Rhombohedral', 'a': 4.51, 'alpha': 57.60             },
    'Te':{'magmom': 2,  'symmetry': 'hcp', 'c/a': 1.330, 'a': 4.45                        },
    'I' :{'magmom': 1,  'symmetry': 'Orthorhombic', 'c/a': 1.347, 'a': 7.27, 'b/a': 0.659 },
    'Xe':{'magmom': 0,  'symmetry': 'fcc', 'a': 6.20                                      },
    'Cs':{'magmom': 1,  'symmetry': 'BCC', 'a': 6.05                                      },
    'Ba':{'magmom': 2,  'symmetry': 'BCC', 'a': 5.02                                      },
    'La':{'magmom': 3,  'symmetry': 'hcp', 'c/a': 1.619, 'a': 3.75                        },
    'Ce':{'magmom': 4,  'symmetry': 'fcc', 'a': 5.16                                      },
    'Pr':{'magmom': 5,  'symmetry': 'hcp', 'c/a': 1.614, 'a': 3.67                        },
    'Nd':{'magmom': 6,  'symmetry': 'hcp', 'c/a': 1.614, 'a': 3.66                        },
    'Pm':{'magmom': 6,                                                                    },
    'Sm':{'magmom': 6,  'symmetry': 'Rhombohedral', 'a': 9.00, 'alpha': 23.13             },
    'Eu':{'magmom': 6,  'symmetry': 'BCC', 'a': 4.61                                      },
    'Gd':{'magmom': 6,  'symmetry': 'hcp', 'c/a': 1.588, 'a': 3.64                        },
    'Th':{'magmom': 6,  'symmetry': 'hcp', 'c/a': 1.581, 'a': 3.60                        },
    'Dy':{'magmom': 6,  'symmetry': 'hcp', 'c/a': 1.573, 'a': 3.59                        },
    'Ho':{'magmom': 6,  'symmetry': 'hcp', 'c/a': 1.570, 'a': 3.58                        },
    'Er':{'magmom': 6,  'symmetry': 'hcp', 'c/a': 1.570, 'a': 3.56                        },
    'Tm':{'magmom': 6,  'symmetry': 'hcp', 'c/a': 1.570, 'a': 3.54                        },
    'Yb':{'magmom': 6,  'symmetry': 'fcc', 'a': 5.49                                      },
    'Lu':{'magmom': 6,  'symmetry': 'hcp', 'c/a': 1.585, 'a': 3.51                        },
    'Hf':{'magmom': 4,  'symmetry': 'hcp', 'c/a': 1.582, 'a': 3.20                        },
    'Ta':{'magmom': 5,  'symmetry': 'BCC', 'a': 3.31                                      },
    'W' :{'magmom': 6,  'symmetry': 'BCC', 'a': 3.16                                      },
    'Re':{'magmom': 5,  'symmetry': 'hcp', 'c/a': 1.615, 'a': 2.76                        },
    'Os':{'magmom': 4,  'symmetry': 'hcp', 'c/a': 1.579, 'a': 2.74                        },
    'Ir':{'magmom': 3,  'symmetry': 'fcc', 'a': 3.84                                      },
    'Pt':{'magmom': 2,  'symmetry': 'fcc', 'a': 3.92                                      },
    'Au':{'magmom': 1,  'symmetry': 'fcc', 'a': 4.08                                      },
    'Hg':{'magmom': 0,  'symmetry': 'Rhombohedral', 'a': 2.99, 'alpha': 70.45             },
    'Tl':{'magmom': 3,  'symmetry': 'hcp', 'c/a': 1.599, 'a': 3.46                        },
    'Pb':{'magmom': 4,  'symmetry': 'fcc', 'a': 4.95                                      },
    'Bi':{'magmom': 3,  'symmetry': 'Rhombohedral', 'a': 4.75, 'alpha': 57.14             },
    'Po':{'magmom': 2,  'symmetry': 'SC', 'a': 3.35                                       },
    'At':{'magmom': 1,                                                                    },
    'Rn':{'magmom': 0,                                                                    },
    'Fr':{'magmom': 1,                                                                    },
    'Ra':{'magmom': 2,                                                                    },
    'Ac':{'magmom': 3,  'symmetry': 'fcc', 'a': 5.31                                      },
    'Th':{'magmom': 4,  'symmetry': 'fcc', 'a': 5.08                                      },
    'Pa':{'magmom': 5,  'symmetry': 'Tetragonal', 'c/a': 0.825, 'a': 3.92                 },
    'U' :{'magmom': 6,  'symmetry': 'Orthorhombic', 'c/a': 2.056, 'a': 2.85, 'b/a': 1.736 },
    'Np':{'magmom': 6,  'symmetry': 'Orthorhombic', 'c/a': 1.411, 'a': 4.72, 'b/a': 1.035 },
    'Pu':{'magmom': 6,  'symmetry': 'Monoclinic'                                          },
    'Am':{'magmom': 6,                                                                    },
    'Cm':{'magmom': 6,                                                                    },
    'Bk':{'magmom': 6,                                                                    },
    'Cf':{'magmom': 6,                                                                    },
    'Es':{'magmom': 6,                                                                    },
    'Fm':{'magmom': 6,                                                                    },
    'Md':{'magmom': 6,                                                                    },
    'No':{'magmom': 6,                                                                    },
    'Lr':{'magmom': 6,                                                                    }}
    

# This data is from Ashcroft and Mermin. --> replaced by hash
reference_states = [\
    None, #X
    {'symmetry': 'diatom', 'd': 0.74}, #H
    {'symmetry': 'atom'}, #He
    {'symmetry': 'BCC', 'a': 3.49}, #Li
    {'symmetry': 'hcp', 'c/a': 1.567, 'a': 2.29}, #Be
    {'symmetry': 'Tetragonal', 'c/a': 0.576, 'a': 8.73}, #B
    {'symmetry': 'Diamond', 'a': 3.57},#C
    {'symmetry': 'diatom', 'd': 1.10},#N
    {'symmetry': 'diatom', 'd': 1.21},#O
    {'symmetry': 'diatom', 'd': 1.42},#F
    {'symmetry': 'fcc', 'a': 4.43},#Ne
    {'symmetry': 'BCC', 'a': 4.23},#Na
    {'symmetry': 'hcp', 'c/a': 1.624, 'a': 3.21},#Mg
    {'symmetry': 'fcc', 'a': 4.05},#Al
    {'symmetry': 'Diamond', 'a': 5.43},#Si
    {'symmetry': 'Cubic', 'a': 7.17},#P
    {'symmetry': 'Orthorhombic', 'c/a': 2.339, 'a': 10.47,'b/a': 1.229},#S
    {'symmetry': 'Orthorhombic', 'c/a': 1.324, 'a': 6.24, 'b/a': 0.718},#Cl
    {'symmetry': 'fcc', 'a': 5.26},#Ar
    {'symmetry': 'BCC', 'a': 5.23},#K
    {'symmetry': 'fcc', 'a': 5.58},#Ca
    {'symmetry': 'hcp', 'c/a': 1.594, 'a': 3.31},#Sc
    {'symmetry': 'hcp', 'c/a': 1.588, 'a': 2.95},#Ti
    {'symmetry': 'BCC', 'a': 3.02},#V
    {'symmetry': 'BCC', 'a': 2.88},#Cr
    {'symmetry': 'Cubic', 'a': 8.89},#Mn
    {'symmetry': 'BCC', 'a': 2.87},#Fe
    {'symmetry': 'hcp', 'c/a': 1.622, 'a': 2.51},#Co
    {'symmetry': 'fcc', 'a': 3.52},#Ni
    {'symmetry': 'fcc', 'a': 3.61},#Cu
    {'symmetry': 'hcp', 'c/a': 1.856, 'a': 2.66},#Zn
    {'symmetry': 'Orthorhombic', 'c/a': 1.695, 'a': 4.51, 'b/a': 1.001},#Ga
    {'symmetry': 'Diamond', 'a': 5.66},#Ge
    {'symmetry': 'Rhombohedral', 'a': 4.13, 'alpha': 54.10},#As
    {'symmetry': 'hcp', 'c/a': 1.136, 'a': 4.36},#Se
    {'symmetry': 'Orthorhombic', 'c/a': 1.307, 'a': 6.67, 'b/a': 0.672},#Br
    {'symmetry': 'fcc', 'a': 5.72},#Kr
    {'symmetry': 'BCC', 'a': 5.59},#Rb
    {'symmetry': 'fcc', 'a': 6.08},#Sr
    {'symmetry': 'hcp', 'c/a': 1.571, 'a': 3.65},#Y
    {'symmetry': 'hcp', 'c/a': 1.593, 'a': 3.23},#Zr
    {'symmetry': 'BCC', 'a': 3.30},#Nb
    {'symmetry': 'BCC', 'a': 3.15},#Mo
    {'symmetry': 'hcp', 'c/a': 1.604, 'a': 2.74},#Tc
    {'symmetry': 'hcp', 'c/a': 1.584, 'a': 2.70},#Ru
    {'symmetry': 'fcc', 'a': 3.80},#Rh
    {'symmetry': 'fcc', 'a': 3.89},#Pd
    {'symmetry': 'fcc', 'a': 4.09},#Ag
    {'symmetry': 'hcp', 'c/a': 1.886, 'a': 2.98},#Cd
    {'symmetry': 'Tetragonal', 'c/a': 1.076, 'a': 4.59},#In
    {'symmetry': 'Tetragonal', 'c/a': 0.546, 'a': 5.82},#Sn
    {'symmetry': 'Rhombohedral', 'a': 4.51, 'alpha': 57.60},#Sb
    {'symmetry': 'hcp', 'c/a': 1.330, 'a': 4.45},#Te
    {'symmetry': 'Orthorhombic', 'c/a': 1.347, 'a': 7.27, 'b/a': 0.659},#I
    {'symmetry': 'fcc', 'a': 6.20},#Xe
    {'symmetry': 'BCC', 'a': 6.05},#Cs
    {'symmetry': 'BCC', 'a': 5.02},#Ba
    {'symmetry': 'hcp', 'c/a': 1.619, 'a': 3.75},#La
    {'symmetry': 'fcc', 'a': 5.16},#Ce
    {'symmetry': 'hcp', 'c/a': 1.614, 'a': 3.67},#Pr
    {'symmetry': 'hcp', 'c/a': 1.614, 'a': 3.66},#Nd
    None,#Pm
    {'symmetry': 'Rhombohedral', 'a': 9.00, 'alpha': 23.13},#Sm
    {'symmetry': 'BCC', 'a': 4.61},#Eu
    {'symmetry': 'hcp', 'c/a': 1.588, 'a': 3.64},#Gd
    {'symmetry': 'hcp', 'c/a': 1.581, 'a': 3.60},#Th
    {'symmetry': 'hcp', 'c/a': 1.573, 'a': 3.59},#Dy
    {'symmetry': 'hcp', 'c/a': 1.570, 'a': 3.58},#Ho
    {'symmetry': 'hcp', 'c/a': 1.570, 'a': 3.56},#Er
    {'symmetry': 'hcp', 'c/a': 1.570, 'a': 3.54},#Tm
    {'symmetry': 'fcc', 'a': 5.49},#Yb
    {'symmetry': 'hcp', 'c/a': 1.585, 'a': 3.51},#Lu
    {'symmetry': 'hcp', 'c/a': 1.582, 'a': 3.20},#Hf
    {'symmetry': 'BCC', 'a': 3.31},#Ta
    {'symmetry': 'BCC', 'a': 3.16},#W
    {'symmetry': 'hcp', 'c/a': 1.615, 'a': 2.76},#Re
    {'symmetry': 'hcp', 'c/a': 1.579, 'a': 2.74},#Os
    {'symmetry': 'fcc', 'a': 3.84},#Ir
    {'symmetry': 'fcc', 'a': 3.92},#Pt
    {'symmetry': 'fcc', 'a': 4.08},#Au
    {'symmetry': 'Rhombohedral', 'a': 2.99, 'alpha': 70.45},#Hg
    {'symmetry': 'hcp', 'c/a': 1.599, 'a': 3.46},#Tl
    {'symmetry': 'fcc', 'a': 4.95},#Pb
    {'symmetry': 'Rhombohedral', 'a': 4.75, 'alpha': 57.14},#Bi
    {'symmetry': 'SC', 'a': 3.35},#Po
    None,#At
    None,#Rn
    None,#Fr
    None,#Ra
    {'symmetry': 'fcc', 'a': 5.31},#Ac
    {'symmetry': 'fcc', 'a': 5.08},#Th
    {'symmetry': 'Tetragonal', 'c/a': 0.825, 'a': 3.92},#Pa
    {'symmetry': 'Orthorhombic', 'c/a': 2.056, 'a': 2.85, 'b/a': 1.736},#U
    {'symmetry': 'Orthorhombic', 'c/a': 1.411, 'a': 4.72, 'b/a': 1.035},#Np
    {'symmetry': 'Monoclinic'},#Pu
    None,#Am
    None,#Cm
    None,#Bk
    None,#Cf
    None,#Es
    None,#Fm
    None,#Md
    None,#No
    None]#Lr

atomic_ichar = [\
    0,#X
    
    1,#H
    0,#He

    1,#Li
    2,#Be
    3,#B
    4,#C
    5,#N
    6,#O
    7,#F
    0,#Ne

    1,#Na
    2,#Mg
    3,#Al
    4,#Si
    5,#P
    6,#S
    7,#Cl
    0,#Ar

    1,#K
    2]#Ca


from math import sin, cos, pi   
reference_lattice = {

    'atom':   [[20.0, 0, 0],
               [0, 20.0, 0],
               [0, 0, 20.0]],

    'diatom': [[20.0, 0, 0],
               [0, 20.0, 0],
               [0, 0, 20.0]],

    'bcc':    [[ 0.5,  0.5, -0.5],
               [-0.5,  0.5,  0.5],
               [ 0.5, -0.5,  0.5]],

    'fcc':    [[ 0.5,  0.5,  0.0],
               [ 0.0,  0.5,  0.5],
               [ 0.5,  0.0,  0.5]],

    'hcp':    [[ 1.0,  0.0,  0.0],
               [ cos(pi/3.), sin(pi/3.),  0.0],
               [ 0.0,  0.0,  1.0]]

}

