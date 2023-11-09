import os, shutil
import argparse
import sys
import time
from oorinano.utils.auxil import convert_time2human as time_convert

from oorinano.calculator.siesta import Siesta
from oorinano.simulators.qtnegf import qtNegf
from oorinano.simulators.qtnegf.cellmodeling import model_electrode, model_channel
from oorinano.simulators.qtnegf.qtplot import qtPlot

def runQtNegf(job, ch_struct, ch_size, el_structs, junc_dist, in_yaml, out_yaml, fdf_params, model_path, nnode, nproc):
    '''
    Input Params
        ch_struct       1 fdf or 3 fdf or key-word to generate
        ch_size         channel size in constructing model
        junc_dist       distance between channel and electrode part in consting model
        el_structs      fdf or atom name to get default
        el_size
        fdf_params       ['elec.fdf', 'scatter.fdf'] to update default parameters
    
    Variables:
        sub_dir         sub directory for qt calculation
        dict_elect      {'psf': path2psf, 'fdf': path2fdf}
        dict_channel    {'psf': path2psf, 'fdf': path2fdf}
     
        model_el        model for electrode, both elements should not be None
        model_structure   model for scattering region, both elements should not be None
    '''
    cwd = os.getcwd()
    
    qt_dir = ['elec', 'channel', 'postprocess']
    #qt_dir = ['.', '.', '3postprocess']
    calc = Siesta()
    ### 1 Make models
    model_el_path = model_path + '/elec'
    model_ch_path = model_path + '/channel'
    dict_elec = model_electrode(calc, el_structs, path2model=model_el_path)
    dict_channel = model_channel(calc, ch_struct, ch_size, junc_dist, path2model=model_ch_path)

    ### 2 Jobs: 
    if job == 'model':
        shutil.copy(dict_elec['fdf'][0], cwd)
        shutil.copy(dict_elec['fdf'][1], cwd)
    elif job == 'run' or job == 'params':
        if job == 'run':
            show_params = False
        else:
            show_params = True
        qtNegf(calc, dict_elec, dict_channel, qt_dir, in_yaml, out_yaml,fdf_params, np=nproc, show_params=show_params)
    elif job == 'plot':         
        #qtPlot(calc, qt_dir[2], qt_dir[1])
        qtPlot(calc, 'image', qt_dir[1], out_yaml)
    return 0

def main():
    cwd = os.getcwd()
    parser = argparse.ArgumentParser(description="run quantum transport at one time")
    parser.add_argument('-j', '--job', default='run', choices=['run', 'model', 'params', 'plot'], help='run qtnegf, make model, check fdf, postprocess only in case "elec", "channel" exist')
    parser.add_argument('-c', '--channel_struct', nargs='*', default=['grp'], help="fdf files: 1 scatter, 2 elect-region, 3 both, or auto-generation")
    parser.add_argument('-cs', '--channel_size', default=6, type=int, help='size of channel')
    parser.add_argument('-mp', '--model_path', default=f'{cwd}/Models', help='if not making models, provide model file path or models')
    gtransport = parser.add_argument_group(title='args for quantum transport for device')
    gtransport.add_argument('-e', '--elec_struct', nargs='*', default=['Au'], help="[one atom to bring fdf|two fdf for right and left]")
    #gtransport.add_argument('-es', '--elec_size', help="to design the size of electrode")
    gtransport.add_argument('-jd', '--junc_dist', default=1.96, type=float, help='distance between electrode & channel')
    gtransport.add_argument('-inpf', '--inp_file', default='input.yaml', help='fname for qt calculation')
    gtransport.add_argument('-outf', '--out_file', default='output.yaml', help='output of qt calculation and input for postprocess')
    gtransport.add_argument('-p', '--params', nargs='*', help='["param_elec.fdf", "param_scat.fdf"] for electrode, scatt calculation parameter updates')
    gprocess = parser.add_argument_group(title='process related arguments')
    gprocess.add_argument("-x", "--partition",  default='X3', help='partition name such as X1, X2, etc')
    gprocess.add_argument("-np", "--nproc",  type=int, default=1, help='number of nprocess')
    gprocess.add_argument("-n", "--nnode",  type=int, default=1, help='number of Nodes')
    parser.add_argument('-u', '--usage', action='store_true', help='usage for run_qt')
    args = parser.parse_args()

    if args.usage:
        print(f"Usage::\
            \n\tThis is how to run quantum transport\
			\n\tCheck 'readme.txt' for detail\
			\n\tRun:\
            \n\t    Direct run\
            \n\t\tpython run_qtnegf.py -m scatt_struct -ms model_size -e elec_struct -es elec_size -jd junc_dist -n nnodes -np nproc\
            \n\t\te.g.:\
            \n\t\t    python run_qtnegf.py -j run -c grp -cs 6 -e Au -jd 1.9 -np 20\
            \n\t\t    python run_qtnegf.py -np 20\
            \n\t\t    python run_qtnegf.py -j model -c grp -cs 6 -e Au -jd 1.9\
            \n\t\t    python run_qtnegf.py -j params -p elec.fdf\
            \n\t\t    * Before run, '$make clean' to delete output subdirectories\
            \n\t\tOptions:\
            \n\t\t    -j    [run (all calculation)|model (print model at workdir)| params (show fdf parameters)]\
            \n\t\t    -c    keyword such as 'grp' or 1 fdf scattering structure and/or psf files\
            \n\t\t    -cs   size of 'grp' in case -c 'grp'\
            \n\t\t    -e   Au| Au.psf and/or left.fdf right.fdf\
            \n\t\t    -p   input parameters in 'param_elec.fdf', 'param_scat.fdf' in wdir\
            \n\t    Job scheduler: slurm\
            \n\t\tsbatch -J test -p {args.partition} -N {args.nnode} -n {args.nproc} slm_qtnegf.sh     # running in main directory\
            \n\t\tsbatch -J test -p {args.partition} -N {args.nnode} -n {args.nproc} --export=sub=0 slm_qtnegf.sh # making subdirectory\
            \n\t\tAs for customized input models, modify slm_qtnegf.sh\
            \n\t\tOptions for sbatch:\
            \n\t\t    -p    partition\
            \n\t\t    -N    number of nodes\
            \n\t\t    -n    number of total process\
            \n\t\t    * Modify python argument inside 'slm_qtnegf.sh'\
			")
        sys.exit(0)
    
    
    start = int(time.time())

    runQtNegf(args.job, args.channel_struct, args.channel_size, args.elec_struct, args.junc_dist, args.inp_file, args.out_file, args.params, args.model_path, args.nnode, args.nproc)

    end  = int(time.time())
    fname = f"lapsetime{args.partition}np{args.nproc}.dat"
    lapse = end - start
    time_str = time_convert(lapse)
    with open(fname, 'w') as f:
        f.write(f"{time_str}")
    return 0

if __name__ == "__main__":
    main()
