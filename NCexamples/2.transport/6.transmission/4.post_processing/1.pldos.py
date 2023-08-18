from nanocore import io
from nanocore.simulator import siesta as s2
from matplotlib import cm
import numpy as np
import time, shutil,os,glob

cwd = os.getcwd()
os.chdir('../3.scatter+tbtrans/siesta')
files = os.listdir()
for f in files:
    if os.path.splitext(f)[-1] in ['.PDOS', '.EIG'] or f == 'STRUCT.fdf':
        shutil.copy(f, cwd)
os.chdir(cwd)

atoms = io.read_struct('STRUCT.fdf') 
#print atoms
# slice atoms by z coordinates
z_coords = []; indice = []
for atom in atoms:
    if not atom[2] in z_coords: z_coords.append(atom[2])

for z in z_coords:
    temp = []
    for atom in atoms:
        if abs(z-atom[2]) < 0.01: temp.append(atom.get_serial())
    indice.append(temp)
#print indice
simobj = 0.0
#find fermilevel
a=glob.glob('*.EIG')
a = str(a[0])
f=open(a)
list_lines=[]
for line in f.readlines():
    list_lines.append(line)
Fermi = list_lines[0]
Fermi = Fermi.split()
Fermi = float(Fermi[0])
print(Fermi)
# get pdos
Z = []; E = []
for ind in indice:
    E1, dos11, dos12 = s2.get_pdos(simobj, -5, 5, by_atom=1, atom_index=ind, broad= 0.05, npoints = 1500, label = 'scatter')
    E = np.array(E1)
    Z.append(np.array(dos11))

absZ = np.abs(Z).T
print(len(E))
# convert to log10 values + minimum correction to avoid -INF
Z = np.log10(absZ)

# generate meshgrid
X, Y = np.meshgrid(np.array(z_coords), E-Fermi)

# customized figure 
import matplotlib.pyplot as plt
fig1 = plt.figure(figsize=(10,5))
# plt.yticks(np.arange(-4,4,1))
levels = np.linspace(1.01*Z.min(), 0.99*Z.max(), 100)
cmap=plt.cm.get_cmap("jet")
import pylab as plb
cset = plb.contourf(X,Y,Z, levels, cmap=cmap)
plb.colorbar(cset,ticks=[-4,-3,-2,-1,0,1])
#plt.show()
fig1.savefig('pldos.png')

