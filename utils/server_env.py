#### server name
import os
import subprocess

###### YHKIM's lab
### partition and number of process per node
np_Xn = { 'X1': 8, 'X2': 12, 'X3': 20, 'X4':24, 'X5':32, 'X6':32 }

#print(f'nproc is {nXn[4]} in {__name__}')

###### KISTI
host_name = subprocess.check_output('hostname', shell=True)
hostname = host_name.decode()
#print(f"host {hostname}")

if  'login' in hostname:
    host = 'kisti'
    home = "/home01/x2818a02"
else: # 'tgm' in hostname:
    host = 'cluster'
    home = "/home/joonho"

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


