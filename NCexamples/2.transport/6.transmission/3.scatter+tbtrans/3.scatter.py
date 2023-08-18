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

voltage = option['Voltage'].split()[0]
if voltage not in os.listdir():
    os.mkdir(voltage)
os.chdir(voltage)
cwd = os.getcwd()

shutil.copytree('../input', 'TSHS')
os.chdir('TSHS')
sim.set_option('TS.Elecs.Neglect.Principal', True)
sim.set_option('kgrid_Monkhorst_Pack', True, (4,1,1))
sim.set_option('PDOS.kgrid_Monkhorst_Pack', True, (12,3,3))
sim.set_option('ProjectedDensityOfStates', True, "-5.00  5.00  0.05  1000 eV\n")
sim.set_option('LocalDensityOfStates', True, "-0.1 0.1 eV \n")
sim.run("scatter", 12, **option)
option['scatter'] = os.getcwd()
os.chdir(cwd)

shutil.copytree('../input', 'TBTrans')
os.chdir('TBTrans')
sim.read_all_fdf()
sim.set_option('kgrid_Monkhorst_Pack', True, (12,1,1))
sim.set_option('TS.Elecs.Neglect.Principal', True)
sim.run("tbtrans", 12, **option)
