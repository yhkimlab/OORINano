from nanocore import siesta as s2
import shutil
import os
import yaml
import argparse

def run_electrode(nproc):
    cwd = os.getcwd()
    #print(f"in 1exec.py cwd {cwd}")
    flist = os.listdir()
    sim = s2.Siesta()
    sim.set_mode('elec')
    for folder in flist:
        if os.path.isdir(folder):
            os.chdir(folder)
            if not os.path.isdir('OUT'):
                shutil.copytree('input', 'OUT')
                os.chdir('OUT')
                sim.read_all_fdf()
                sim.run(nproc)
            os.chdir(cwd)
    return 0            

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--nproc", type=int, default=1)
    args = parser.parse_args()

    run_electrode(args.nproc)
    return 0

if __name__ == '__main__':
    main()
