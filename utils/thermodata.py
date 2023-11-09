'''
Data to be improted for catalysis
Units   eV
Future
    DFT functional dependance can be included
'''

### PBE
Etot_H2O    = -14.236           # H2O(g) in VASP
Etot_H2     =  -6.760           # H2(g)  in VASP

### Norskov: (2004, table 3)
zpe_H2O      = 0.560            
zpe_H2       = 0.270             

### NIST: https://janaf.nist.gov/
S_H2O       =  0.001957         # NIST 0.58; it will be 0.67 at P=0.035 bar
S_H2        = 0.001354          # NIST 0.40; Norskov 0.41


