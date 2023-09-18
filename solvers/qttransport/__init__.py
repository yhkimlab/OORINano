from __future__ import print_function
from ...atoms   import *
from ...aux.env import *
import sys, os, shutil, importlib     ### for 1elec directory

from nanocore.models import carbonlab
from nanocore.aux.convert import convert_xyz2abc
from nanocore.ncio import write_struct, read_struct ### for 2model

"""
Transport(elec1, chan, elec2=None)

Class for representing a quantum transport simulation

Parameters
----------
electrode1/2 : AtomsSystem
    an atomic structure whose PBC is imposed along the transport direction
channel : AtomsSystem 
    an atomic structure 

Optional parameters
-------------------

Examples
--------
>>> elec = carbonlab.grp_armchair(4,4)
>>> chan = carbonlab.chain(10, 1.285, 0.1)
>>> trans_obj = Transport(elec, chan)


__slots__ = ['_elec1', '_elec2', '_chan']

"""

elec_atom = ['Au']

def qtmodel(model_str, chsize, el_str, junc_dist):


def get_elec_inifile(calc, files):
    '''
    to get default pathway, calc is passed
    '''
    cwd = os.getcwd()
    modulename = importlib.import_module(calc.__class__.__module__)
    ### .psf file
    if type(files) == str:
        dir_default = modulename.default_location+'/psf'
        fpath = dir_default+'/'+files
        if os.path.isfile(fpath):
            return fpath
        else:
            print(f"Can't find {fpath}")
            sys.exit(101)
    ### electrode 2 fdf files
    else:
        dir_default = modulename.default_location+'/model/elec'
        fpaths=[]
        for f in files:
            fpath = dir_default+'/'+f
            if os.path.isfile(fpath):
                fpaths.append(fpath)
            else:
               print(f"Can't find {fpath}")
               sys.exit(102)
        return fpaths
                


def qtTransport(job, calc, junc_dist=None, np=1, det_k=0, max_nk=10, init_nk=1, step_nk=1, opts=None):

    ### 1. Electrode calculation
    calc_electrode(job, calc, np)
    ### 2. Make model for scattering region
    make_scattering_model(job, calc, junc_dist)

    return 0
    
def make_scattering_model(job, calc, junc_dist ):
    cwd = os.getcwd()
    modulename = importlib.import_module(calc.__class__.__module__)
    dir_default = modulename.default_location+'/model/scatter'
    ###### get .psf files
    path_pp=[]
    for at in calc.model_scatter[0]:
        if os.path.splitext(at)[-1] == '.psf':
            path_psf = cwd+'/'+at
        else:
            f_psf = at+'.psf'
            path_psf = modulename.default_location+f'/psf/{f_psf}'
        path_pp.append(path_psf)
    
    ###### get fdf file
    ### calc.model_scatter[1]: [channel, left-elect, right-elect]
    ### if only channel given, get default electrode part
    if len(calc.model_scatter[1]) <= 1:
        path_scatter_left = dir_default+'/elecL.fdf'
        path_scatter_rigt = dir_default+'/elecR.fdf'
    ### if electrode part exists
    else:
        path_scatter_left = cwd + f'/calc.model_scatter[1][1]'
        path_scatter_rigt = cwd + f'/calc.model_scatter[1][2]'

    elecL = read_struct(path_scatter_left)
    elecR = read_struct(path_scatter_rigt)

    #for i in range(min_unit, max_unit, 1):
    cnt = carbonlab.grp_rect(2,i)
    cnt.select_all()
    cnt.rotate(90, axis_dir=(1,0,0), with_cell=True)
    cnt.translate(-1.42, 0, 0)
    cnt.wrap_positions()
    cnt.select_z(-0.1, 1)
    cnt.delete()

    cnt_cell = cnt.get_cell()
    elecL_cell = elecL.get_cell()
    elecL_cell_abc = convert_xyz2abc(elecL_cell[0], elecL_cell[1], elecL_cell[2])
    cnt_cell_abc = convert_xyz2abc(cnt_cell[0], cnt_cell[1], cnt_cell[2])
    ratio = elecL_cell_abc[0] / cnt_cell_abc[0]
    cnt2 = cnt.adjust_cell_size(ratio, 7)

    cnt2.select_all()
    cnt2.translate(0,5,junc_dist+elecL.get_zmax()-cnt2.get_zmin())
    new = elecL + cnt2
    elecR.select_all()
    elecR.translate(0,0,junc_dist+new.get_zmax()-elecR.get_zmin())
    new = new + elecR
    elecL_cell[2][2] = new.get_zmax() - new.get_zmin()
    new.set_cell(elecL_cell)
    new.set_vacuum(2.306)

    new_cell = new.get_cell()
    new_zavg = (new.get_zmin()+new.get_zmax())/2
    cell_avg = convert_xyz2abc(new_cell[0], new_cell[1], new_cell[2])[2]/2
    new.select_all()
    new.translate(0,0,cell_avg - new_zavg)
    new.select_all()
    new.sort(option = 'z'); new.set_serials(1)
    write_struct(new, fname = f"cnt_{i}.fdf")
    # write_poscar(new, f"cnt_{i}.poscar")


    return 0
    
def calc_electrode(job, calc, np):
    calc.set_mode('elec')
    
    cwd = os.getcwd()
    

    ###### 1 Make electrode model 
    ### if psf, else get default M.psf
    if os.path.splitext(calc.model_el[0])[-1] == '.psf':
        path_psf = cwd+'/'+calc.model_el[0]
    elif calc.model_el[0] in elec_atom:
        f_psf = calc.model_el[0]+'.psf'
        path_psf = get_elec_inifile(calc, f_psf)
    ### if None, get default 
    if calc.model_el[1] is None:
        f_fdf_1 = calc.model_el[0]+'_STRUCT_left.fdf'
        f_fdf_2 = calc.model_el[0]+'_STRUCT_rigt.fdf'
        path_fdf_left, path_fdf_rigt = get_elec_inifile(calc, [f_fdf_1, f_fdf_2])
    else:
        path_fdf_left = cwd+'/'+calc.model_el[1][0]
        path_fdf_rigt = cwd+'/'+calc.model_el[1][1]
        
    if job == 'model':
        print(f"copy {path_fdf_left} to {cwd}")
        shutil.copy(path_fdf_left, cwd)
        shutil.copy(path_fdf_rigt, cwd)
        return 0
    
    ### Make main subdirectory for electrode calculation

    if not os.path.exists(calc.sub_dir[0]):
        os.mkdir(calc.sub_dir[0])
        os.chdir(calc.sub_dir[0])
        ### two calculations for left and right electrode
        ### make two subdirectories inside 1elec dir
        sub_dirs=[]
        sub_dirs.append(calc.model_el[0]+'_left')
        sub_dirs.append(calc.model_el[0]+'_rigt')
        fdf_files = [path_fdf_left, path_fdf_rigt]
        
        for subdir, fstruct in zip(sub_dirs, fdf_files):
            os.mkdir(subdir)
            os.chdir(subdir)
            os.mkdir('Input')
            shutil.copy(path_psf, 'Input')
            shutil.copy(fstruct, 'Input/STRUCT.fdf')
            if job == 'run':
                shutil.copytree('Input', 'Run')
                os.chdir('Run')
                calc.read_all_fdf()
                calc.run_qt(np)
                print(f"Job completed in {subdir}")
                os.chdir('..')
            os.chdir('..')
    os.chdir(cwd)
    return 0

        
        
    '''
        = get_elect_inifile(cwd)
        os.mkdir('input')



    ### check input file
    # 1. set electrode and channel
    set_electrode1(elec1)
    if not elec2: elec2 = elec1.copy()
    set_electrode2(elec2)
    set_channel(chan)

    # 2. check: cell shape & size
    possible = self.is_possible()
    #if not possible: raise ValueError, "Please check your initial structures."

    # 3. determine number of kpoint
    if det_k: 
        kpts, energies = self.determine_nkpt(self._elec1,
                                             max_nk, init_nk, step_nk, opts)
        i = 0
        for kpt in kpts:
            print (kpt, energies[i])
            i += 1

    # 4. make junction model
    #self.make_junction(auto_adjust)
    '''

