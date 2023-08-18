import os, sys, shutil
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from nanocore.simulator.siesta import get_hartree_pot_z

cwd = os.getcwd()
os.chdir('../3.scatter+tbtrans/0.0/TSHS')
files = [f for f in os.listdir() if os.path.splitext(f)[-1] == '.VH']
assert len(files) == 1
shutil.copy(files[0], cwd)
os.chdir(cwd)

pos, pot = get_hartree_pot_z(os.path.splitext(files[0])[0])
pos = np.asarray(pos); pos = pos * 0.5291
fig = plt.figure(figsize=[7,4])
fig.suptitle('Potential', fontsize = 30)
ax1 = fig.add_subplot(1,1,1)

ax1.plot(pos, pot,'r', linewidth=4)
ax1.set_xlabel('Z [Ang]',fontsize = 15)
ax1.set_ylabel('$V_{H} [eV]$',fontsize = 15)
ax1.tick_params(axis='both' ,labelsize=15, length = 6, width = 3, direction='in')
plt.setp(ax1.spines.values(), linewidth=4)
[ax1.spines[i].set_linewidth(3) for i in ax1.spines.keys()]
ax1.tick_params(axis='both',which='major', labelsize=15,direction='in', length = 12, width = 3)
ax1.tick_params(axis='both',which='minor', labelsize=15,direction='in', length = 6, width = 3)
ax1.ticklabel_format(useOffset=False)
fig.tight_layout()

plt.savefig('potential.jpg')