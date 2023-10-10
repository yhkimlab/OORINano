import os, shutil
import argparse
import sys
import time
from nanocore.aux import convert_time2human as time_convert

from nanocore.simulators.siesta import Siesta
from nanocore.solvers.qtnegf import qtNegf
from nanocore.solvers.qtnegf.cellmodeling import model_electrode, model_channel
from nanocore.solvers.qtnegf.qtplot import qtPlot

def runQtNegf(job, ch_struct, ch_size, el_structs, el_size, junc_dist, in_yaml, out_yaml, fdf_params, nnode, nproc):

    '''
    Input Params
        ch_struct       1 fdf or 3 fdf or key-word to generate
        ch_size         channel size in constructing model
        junc_dist       distance between channel and electrode part in consting model
        el_structs      fdf or atom name to get default
        el_size
        fdf_params       ['el[ec].fdf', 'scat[tering].fdf'] to update default parameters
    
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
    dict_elec = model_electrode(calc, el_structs, el_size)
    dict_channel = model_channel(calc, ch_struct, ch_size, junc_dist)

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
    elif job == 'paramsss':
        print("=== Electrode Calculation ===")
        calc.set_mode('elec')
        if 'elec.fdf' in fdf_params:
            calc.add_fdf('elec.fdf')
        files = ['BASIS.fdf', 'KPT.fdf', 'RUN.fdf']
        for f in files:
            calc.print_fdf(f)
        print("\n=== Scattering Calculation ===")
        calc.set_clean()
        calc.set_mode('scatter')
        files.append('TS.fdf')
        if 'scatter.fdf' in fdf_params:
            calc.add_fdf('scatter.fdf')
        for f in files:
            calc.print_fdf(f)
    elif job == 'plot':         
        #qtPlot(calc, qt_dir[2], qt_dir[1])
        qtPlot(calc, cwd, qt_dir[1], out_yaml)
    return 0

def main():
    parser = argparse.ArgumentParser(description="run quantum transport at one time")
    parser.add_argument('-j', '--job', default='run', choices=['run', 'model', 'params'], help='run qt, make model, check fdf')
    parser.add_argument('-c', '--channel_struct', nargs='*', default=['grp'], help="fdf files: 1 scatter, 2 elect-region, 3 both, or auto-generation")
    parser.add_argument('-cs', '--channel_size', default=6, type=int, help='size of channel')
    gtransport = parser.add_argument_group(title='args for quantum transport for device')
    gtransport.add_argument('-e', '--elec_struct', nargs='*', default=['Au'], help="[one atom to bring fdf|two fdf for right and left]")
    gtransport.add_argument('-es', '--elec_size', help="to design the size of electrode")
    gtransport.add_argument('-jd', '--junc_dist', default=1.96, type=float, help='distance between electrode & channel')
    gtransport.add_argument('-inpf', '--inp_file', default='input.yaml', help='fname for qt calculation')
    gtransport.add_argument('-outf', '--out_file', default='output.yaml', help='output of qt calculation and input for postprocess')
    gtransport.add_argument('-p', '--params', nargs='*', help='["elec.fdf", "scatt.fdf"] for electrode, scatt calculation parameter updates')
    gprocess = parser.add_argument_group(title='process related arguments')
    gprocess.add_argument("-np", "--nproc",  type=int, default=1, help='number of nprocess')
    gprocess.add_argument("-n", "--nnode",  type=int, default=1, help='number of Nodes')
    parser.add_argument('-u', '--usage', action='store_true', help='usage for run_qt')
    args = parser.parse_args()

    if args.usage:
        print(f"Usage::\
            \n\tThis is how to run quantum transport\
			\n\tCheck 'readme.txt' for detail\
			\n\tRun:\
            \n\t    python run_qtnegf.py -m scatt_struct -ms model_size -e elec_struct -es elec_size -jd junc_dist -n nnodes -np nproc\
            \n\t    e.g.:\
            \n\t\t python run_qtnegf.py -j run -c grp -cs 6 -e Au -jd 1.9 -np 20\
            \n\t\t python run_qtnegf.py -j model -c grp -cs 6 -e Au -jd 1.9\
            \n\t\t python run_qtnegf.py -j params -p elec.fdf\
            \n\t    Options:\
            \n\t\t-j    [run (all calculation)|model (print model at workdir)| params (show fdf parameters)]\
            \n\t\t-c   atoms in scattering model\
            \n\t\t-cs   fdf scattering structure: 1 for channel, 2 for parts of electrode\
            \n\t\t-e   one atom for electrode\
            \n\t\t-es   fdf 2 electrode structure for left & right\
            \n\t\t-p   input parameters in 'elec.fdf', 'scatter.fdf' in wdir\
			")
        sys.exit(0)
    
    
    start = int(time.time())

    runQtNegf(args.job, args.channel_struct, args.channel_size, args.elec_struct, args.elec_size, args.junc_dist, args.inp_file, args.out_file, args.params, args.nnode, args.nproc)

    end  = int(time.time())
    fname = f"lapsetime{args.nproc}.dat"
    lapse = end - start
    time_str = time_convert(lapse)
    with open(fname, 'w') as f:
        f.write(f"{time_str}")
    return 0

if __name__ == "__main__":
    main()
