from nanocore.simulator import siesta as s2
import shutil
import os
import yaml
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--node", dest="node", action="store", type=int, default=1)
args = parser.parse_args()

cwd = os.getcwd()
flist = os.listdir()
sim = s2.Siesta()
sim.set_mode('elec')
for folder in flist:
    if os.path.isdir(folder):
        os.chdir(folder)
        shutil.copytree('input', 'OUT')
        os.chdir('OUT')
        sim.read_all_fdf()
        sim.run(args.node)
        os.chdir(cwd)
