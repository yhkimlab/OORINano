from nanocore.simulator import siesta as s2
import shutil
import os, sys
import yaml
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", dest="input", action="store", type=str, default="input.yaml")
parser.add_argument("-o", "--out", dest="output", action="store", type=str, default="output.yaml")
parser.add_argument("-m", "--model", dest="model", action="store", type=str, default="cnt_6.fdf")
parser.add_argument("-n", "--node", dest="node", action="store", type=int, default=1)
args = parser.parse_args()

model = args.model

cwd = os.getcwd()
shutil.copy(f'../2.model/{model}', 'input/STRUCT.fdf')
with open(args.input) as f:
    option = yaml.safe_load(f)

os.chdir('input')
sim = s2.Siesta()
sim.read_all_fdf()
os.chdir('..')

voltage = option['Voltage'].split()[0]
if voltage not in os.listdir():
    os.mkdir(voltage)
os.chdir(voltage)
cwdv = os.getcwd()

shutil.copytree('../input', 'TSHS')
os.chdir('TSHS')
sim.set_mode('scatter')
sim.set_option('TS.Elecs.Neglect.Principal', True)
sim.run(args.node, **option)
option['scatter'] = os.getcwd()
os.chdir(cwdv)

shutil.copytree('../input', 'TBTrans')
os.chdir('TBTrans')
sim.set_mode('tbtrans')
sim.read_all_fdf()
sim.set_option('kgrid_Monkhorst_Pack', True, (12,1,1))
sim.set_option('TS.Elecs.Neglect.Principal', True)
sim.run(args.node, **option)
option['tbtrans'] = os.getcwd()
os.chdir(cwd)
with open(args.output, 'w') as f:
    yaml.safe_dump(option, f)
