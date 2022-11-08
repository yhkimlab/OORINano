from nanocore import io
from nanocore import catalysis     

at = io.read_poscar('POSCAR')

catalysis.runORR(at, nproc=24, npar=4, mode='opt', kpoints=[4,4,1], vib=1, label='test')

### belows for explicit plot
### if each image is calculated in the separate directory, the calculation can be skipped if already calculation was done
#TE, ZPE, TS = catalysis.runORR(at, nproc=24, npar=4, mode='opt', kpoints=[4,4,1], vib=1, label='test')
#print(TE)
#print(ZPE)
#print(TS)

#G_ORR      = catalysis.gibbs_ORR_4e_acid(TE=TE, pH=0)
#print(G_ORR)
#G_ORR_vib  = catalysis.gibbs_ORR_4e_acid(TE=TE, ZPE=ZPE, TS=TS, pH=0)
#print(G_ORR_vib)
#catalysis.plot_ORR_4e_acid(G_ORR_vib, U=0.7, legend=['U=1.23V', 'U=0.70V', 'U=0.00V'])

