from nanocore.simulator import siesta as s2
import shutil
import os, sys
import yaml

model = "cnt_6.fdf"

shutil.copy(f'../2.model/{model}', 'input/STRUCT.fdf')
with open("input.yaml") as f:
    option = yaml.safe_load(f)

os.chdir('input')
sim = s2.Siesta()
sim.read_all_fdf()
os.chdir('..')


if 'siesta' not in os.listdir():
    shutil.copytree('input', 'siesta')
os.chdir('siesta')
cwd = os.getcwd()
# sim.set_option('TS.Elecs.Neglect.Principal', True)
sim.set_option('kgrid_Monkhorst_Pack', True, (4,1,1))
sim.set_option('PDOS.kgrid_Monkhorst_Pack', True, (12,3,3))
sim.set_option('ProjectedDensityOfStates', True, "-12.00  2.00  0.05  1000 eV\n")
sim.set_option('LocalDensityOfStates', True, "-4.268524 -4.068524 eV \n")
sim.run("siesta", 12, **option)
