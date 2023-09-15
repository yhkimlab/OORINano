import subprocess
import os, shutil
import argparse
import sys
import time
from nanocore.aux.convert import convert_to_preferred_format as time_convert

def runQuantumTransport(nproc):
    cwd = os.getcwd()
    #print(f"in running script, cwd {cwd}")
    os.chdir('1.elec')
    cmd = f'python ../1elec_calc.py -np {nproc}'
    result = subprocess.run(cmd, shell=True, check=True)
    print("Complete electrode calculation.")
    os.chdir(cwd)
    os.chdir('2.model')
    cmd = f'python ../2generate_model.py'
    result = subprocess.run(cmd, shell=True, check=True)
    os.chdir(cwd)
    os.chdir('3.scatter_tbtrans')
    cmd = f'python ../3scatter.py -np {nproc}'
    result = subprocess.run(cmd, shell=True, check=True)
    print("Complete quantum transport calculation")
    shutil.copy("output.yaml", '../4.post_processing')
    os.chdir(cwd)
    os.chdir('4.post_processing')
    cmd = f'python ../4post_process.py'
    result = subprocess.run(cmd, shell=True, check=True)
    os.chdir(cwd)
    return 0

def main():
    parser = argparse.ArgumentParser(description="run quantum transport at one time")
    parser.add_argument("-np", "--nproc",  type=int, default=1, help='number of nprocess')
    parser.add_argument("-n", "--nnode",  type=int, default=1, help='number of Nodes')
    parser.add_argument('-u', '--usage', action='store_true', help='usage for run_qt')
    args = parser.parse_args()

    if args.usage:
        print(f"Usage::\
            \n\tThis is how to run quantum transport\
			\n\tCheck 'readme.txt' for detail\
			\n\t1. direct run\
			\n\t    python run_qt.py -np nproc\
			")
        sys.exit(0)

    start = int(time.time())
    runQuantumTransport(args.nproc)
    end  = int(time.time())
    fname = f"lapsetime{args.nproc}.dat"
    lapse = end - start
    time_str = time_convert(lapse)
    with open(fname, 'w') as f:
        f.write(f"{time_str}")
    return 0

if __name__ == "__main__":
    main()
