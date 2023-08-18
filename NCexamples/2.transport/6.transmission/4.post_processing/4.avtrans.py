from nanocore import io, vis
from nanocore.simulator import siesta as s2
import matplotlib.pyplot as plt
import numpy as np
import time, shutil,os,glob

cwd = os.getcwd()
os.chdir('../3.scatter+tbtrans/0.5/TBTrans')
files = os.listdir()
for f in files:
    if os.path.splitext(f)[-1] in ['.AVTRANS_Left-Right']:
        file = f
        shutil.copy(f, cwd)
os.chdir(cwd)

energy, trans = s2.get_transmission(file)
# trans = np.log10(trans)

fig = plt.figure(figsize=[7,4])
# fig.suptitle('Potential', fontsize = 30)
ax1 = fig.add_subplot(1,1,1)
ax1.plot(energy, trans,'r', linewidth=4)
ax1.set_yscale('log')
ax1.set_xlabel('$E-E_F [eV]$',fontsize = 15)
ax1.set_ylabel('Transmission',fontsize = 15)
ax1.set_ylim(0.01, 1)
ax1.tick_params(axis='both' ,labelsize=15, length = 6, width = 3, direction='in')
plt.setp(ax1.spines.values(), linewidth=4)
[ax1.spines[i].set_linewidth(3) for i in ax1.spines.keys()]
ax1.tick_params(axis='both',which='major', labelsize=15,direction='in', length = 12, width = 3)
ax1.tick_params(axis='both',which='minor', labelsize=15,direction='in', length = 6, width = 3)
fig.tight_layout()

plt.savefig('trans.jpg')
