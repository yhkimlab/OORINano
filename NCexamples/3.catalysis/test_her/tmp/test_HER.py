from NanoCore import surflab
from NanoCore import catalysis

at = surflab.fccsurfaces('Pt', '111', (3,3,3), vac=15)

catalysis.runHER(at, nproc=20, npar=4, kpoints=[4,4,1],ediffg = -0.04)




