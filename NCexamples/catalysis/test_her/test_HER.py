from nanocore import surflab
from nanocore import catalysis
from pyslurm import job

### make slab model
at = surflab.fccsurfaces('Pt', '111', (3,3,3), vac=15)

### run HER: plot in runHER
catalysis.runHER(at, nproc=48, npar=8, kpoints=[4,4,1],ediffg = -0.04)




