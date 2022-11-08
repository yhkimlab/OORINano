from NanoCore import io
from NanoCore import catalysis     

at = io.read_poscar('CONTCAR_Pt-SAC')

catalysis.runORR(at, nproc=36, npar=6, mode='opt', kpoints=[2,2,1], vib=1, label='test')


