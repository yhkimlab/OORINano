from nanocore.simulator import siesta as s2
import shutil
import os
import yaml

with open("input.yaml") as f:
    option = yaml.safe_load(f)
cwd = os.getcwd()
flist = os.listdir()
sim = s2.Siesta()
for folder in flist:
    if os.path.isdir(folder):
        os.chdir(folder)
        shutil.copytree('input', 'OUT')
        os.chdir('OUT')
        sim.read_all_fdf()
        sim.set_option("kgrid_Monkhorst_Pack", True, (4,1,50))
        sim.run("elec", 12, **option)
        os.chdir(cwd)
