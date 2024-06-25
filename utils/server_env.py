#### server name
import os
import subprocess

###### KHKIM's lab
### partition and number of process per node
nXn = { 1: 8, 2: 12, 3: 20, 4:24, 5:32, 6:32 }
#print(f'nproc is {nXn[4]} in {__name__}')

###### KISTI
_HOSTname = subprocess.check_output('hostname', shell=True)
if _HOSTname == 'login02':
    host = 'kisti'
else:
    host = _HOSTname

if host == 'kisti':    
    home = "/home01/x2462a02"

'''
PPNs={'chi':1, 'login':4, 'iron':''}
PBS_QChem_scripts={'chi':'None', 'login':'', 'iron':''}
PBS_VASP_scripts={'chi':'None', 'login':'', 'iron':''}

class Machines:
    """ info for KAIST server """
    def __init__(self, ppn, qc_pbs, vasp_pbs):
        self.ppn = ppn
        self.qc_pbs = qc_pbs
        self.vasp_pbs = vasp_pbs

_HOST = _HOSTname.split()[0].decode('utf-8')
_USER = os.getenv('USER')

home = "/home/"+ _USER
python      = home + '/anaconda3/bin/python'
PBS_HOME    = home + '/sandbox/pypbs/'    

pbs_vname   = PBS_VASP_scripts[_HOST]
pbs_qname   = PBS_QChem_scripts[_HOST]
host_machine=Machines(PPNs[_HOST], PBS_HOME + pbs_qname, PBS_HOME + pbs_vname)
'''


