import subprocess
import os, shutil
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--node", dest="node", action="store", type=int, default=1)
args = parser.parse_args()

cwd = os.getcwd()
os.chdir('1.elec')
cmd = f'python 1.elec_calc.py -n {args.node}'
result = subprocess.run(cmd, shell=True, check=True)
os.chdir('../2.model')
cmd = f'python 2.generate_model.py'
result = subprocess.run(cmd, shell=True, check=True)
os.chdir('../3.scatter+tbtrans')
cmd = f'python 3.scatter.py -n {args.node}'
result = subprocess.run(cmd, shell=True, check=True)
shutil.copy("output.yaml", '../4.post_processing')
os.chdir('../4.post_processing')
cmd = f'python 4.post_process.py'
result = subprocess.run(cmd, shell=True, check=True)