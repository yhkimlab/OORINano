import subprocess
import os, shutil
import argparse
import sys
import time
from nanocore.aux import convert_time2human as time_convert

### OOP and modularize
from nanocore.solvers.qttransport import qtTransport, qtmodel
from nanocore.solvers.qttransport import qtPlot
from nanocore.simulator.siesta import Siesta

def runQuantumTransport(job, model_structs, ch_size, el_structs, el_size, junc_dist, in_yaml, out_yaml, nnode, nproc):

    '''
    model_structs   1 fdf or 3 fdf or key-word to generate
    ch_size         channel size in constructing model
    junc_dist       distance between channel and electrode part in consting model
    el_structs      fdf or atom name to get default
          
    Siesta: 
        model_el        model for electrode, both elements should not be None
        model_structure   model for scattering region, both elements should not be None
        sub_dir         sub directory for qt calculation
    '''
    cwd = os.getcwd()

    qt_dir = ['1elec', '2channel', '3postprocess']
    #qt_dir = ['.', '.', '3postprocess']
    calc = Siesta()
    ### 1 Make models
    dict_elec, dict_model = qtmodel(calc, model_structs, ch_size, el_structs, el_size, junc_dist)
    if job == 'model':
        shutil.copy(dict_elec['fdf'][0], cwd)
        shutil.copy(dict_elec['fdf'][1], cwd)
    elif job == 'run':
        qtTransport(calc, dict_elec, dict_model, qt_dir, in_yaml, out_yaml, np=nproc)
    elif job == 'plot':         
        qtPlot(calc, qt_dir[2], qt_dir[1])
    #elif job == 'fdf':
    #    show_ini
    '''
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
    '''
    return 0

def main():
    parser = argparse.ArgumentParser(description="run quantum transport at one time")
    parser.add_argument('-j', '--job', default='run', choices=['run', 'model', 'fdf'], help='run quantum transport, make model, check fdf')
    parser.add_argument('-m', '--model_struct', nargs='*', default=['grp'], help="model structure or keyword for auto-generation")
    parser.add_argument('-ms', '--model_size', default=6, type=int, help='size of channel')
    gtransport = parser.add_argument_group(title='args for quantum transport for device')
    gtransport.add_argument('-e', '--elec_struct', nargs='*', default=['Au'], help="[one atom to generate structure|two structures file for left(first) and right]")
    gtransport.add_argument('-es', '--elec_size', help="to design the size of electrode")
    gtransport.add_argument('-jd', '--junc_dist', default=1.96, type=float, help='distance between electrode & channel')
    gtransport.add_argument('-inpf', '--inp_file', default='input.yaml', help='fname for qt calculation')
    gtransport.add_argument('-outf', '--out_file', default='output.yaml', help='output of qt calculation and input for postprocess')
    gprocess = parser.add_argument_group(title='process related arguments')
    gprocess.add_argument("-np", "--nproc",  type=int, default=1, help='number of nprocess')
    gprocess.add_argument("-n", "--nnode",  type=int, default=1, help='number of Nodes')
    parser.add_argument('-u', '--usage', action='store_true', help='usage for run_qtransport')
    args = parser.parse_args()

    if args.usage:
        print(f"Usage::\
            \n\tThis is how to run quantum transport\
			\n\tCheck 'readme.txt' for detail\
			\n\tRun:\
            \n\t    python run_qtransport.py -j run -m scatt_struct -ms model_size -e elec_struct -es elec_size -jd junc_dist -n nnodes -np nproc\
            \n\t    python run_qtransport.py -j run -m grp -e Au -jd 1.9 -np 24\
            \n\t    python run_qtransport.py -j model -m grp -e Au -jd 2.5\
            \n\t    Options:\
            \n\t\t-j    [run (calculation)|model (print model at workdir)|fdf (show fdf parameters)]\
            \n\t\t-m    model key word 'grp'\
            \n\t\t-e    one atom species for electrode or two structure files\
            \n\t\t-jd   junction distance\
			")
        sys.exit(0)
    
    
    start = int(time.time())

    runQuantumTransport(args.job, args.model_struct, args.model_size, args.elec_struct, args.elec_size, args.junc_dist, args.inp_file, args.out_file, args.nnode, args.nproc)

    end  = int(time.time())
    fname = f"lapsetime{args.nproc}.dat"
    lapse = end - start
    time_str = time_convert(lapse)
    with open(fname, 'w') as f:
        f.write(f"{time_str}")
    return 0

if __name__ == "__main__":
    main()
