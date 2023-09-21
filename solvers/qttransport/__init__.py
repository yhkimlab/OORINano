from __future__ import print_function
from ...atoms   import *
from ...aux.env import *
import sys, os, re, shutil, importlib     

from nanocore.models import carbonlab
from ...aux   import convert_xyz2abc, check_file ### for model, 1elec

import yaml             ### added for cal scattering
from .qtplot import qtPlot


"""
Transport(elec1, chan, elec2=None)

Class for representing a quantum transport calculation

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
job_complete = "Job completed"

def get_elec_inifile(calc, files):
    '''
    to get default pathway, calc is passed
    '''
    cwd = os.getcwd()
    modulename = importlib.import_module(calc.__class__.__module__)
    
    ### .psf file
    if type(files) == str and re.search('psf', files):
        dir_default = modulename.default_location+'/psf'
        fpath = dir_default+'/'+files
        if os.path.isfile(fpath):
            return fpath
        else:
            print(f"Can't find {fpath}")
            sys.exit(111)
    ### electrode 2 fdf files
    elif type(files) == list:
        dir_default = modulename.default_location+'/model/elec'
        fpaths=[]
        for f in files:
            fpath = dir_default+'/'+f
            if os.path.isfile(fpath):
                fpaths.append(fpath)
            else:
               print(f"Can't find {fpath}")
               sys.exit(112)
        return fpaths
    else:
        print("Input error in ge_elec_inifile")
        sys.exit(103)
        

def qtmodel(calc, model_struct, chsize, el_struct, el_size, junc_dist):
    '''
    calc                to obtain default directory pathway
    To make electrode model
        el_struct       M or M.spf and/or 2 fdf file
    '''
    cwd = os.getcwd()
    
    modulename = importlib.import_module(calc.__class__.__module__)
    
    ###### 1. Make electrode model
        
    fdf = [ f for f in el_struct if 'fdf' in f]
    lfdf = False
    lpsf = False
    ### if Au_STRUCT_left[rigt].fdf exists
    if len(fdf) == 2:
        path_fdf_left = cwd+'/'+fdf[0]
        path_fdf_rigt = cwd+'/'+fdf[1]
        lfdf = True
        el_struct = [ f for f in el_struct if f not in fdf ]
    elif len(fdf) == 1:
        print("Can't run with one-side electrode structure")
        sys.exit(101)

    ### if metal species exists
    if len(el_struct) == 1:
        if os.path.splitext(el_struct[0])[-1] == '.psf':
            path_psf = cwd+'/'+el_struct[0]
            metal = os.path.splitext(el_struct[0])[0]
            lpsf = True
        else:
            metal = el_struct[0]
            f_psf = metal + '.psf'
            
    ### if not metal species
    else:
        metal = modulename.read_struct(path_fdf_left)[0].get_symbol()
        f_psf = metal+'.psf'
    
    ### if not path_psf, get default
    if not lpsf and 'f_psf' in locals():
        path_psf = get_elec_inifile(calc, f_psf)
    else:
        print("Algorithm error")
        sys.exit(102)
    ### if no path_fdf, get default
    if not lfdf and 'metal' in locals():
        f_fdf1 = metal+'_STRUCT_left.fdf'
        f_fdf2 = metal+'_STRUCT_rigt.fdf'
        path_fdf_left, path_fdf_rigt = get_elec_inifile(calc, [f_fdf1, f_fdf2])
    
    dict_elec={'psf': path_psf, 'fdf': [path_fdf_left, path_fdf_rigt]}
    
    
    ######### 2. Make scattering region ################################
    ### 2.1 Generate scattering region

    ###### get fdf file: [channel, left-elect, right-elect]
    ### if only channel given, get default electrode part
    dir_default = modulename.default_location+'/model/scatter'
    
    lfdf = False
    lscatter = False
    lsidepart = False

    fdf = [f for f in model_struct if 'fdf' in f]
    ### No inupt fdf but model name, gnenerate
    #print(f"{fdf}")
    if len(fdf) == 1:
        lfdf = True
        pass
    elif len(fdf) == 0:
        ### if input is one, it is channel name
        if len(model_struct) == 1:
            channel = model_struct[0]
        else:
            print(f"Algorithm error in preparing scattering region")
            sys.exit(110)
    elif len(fdf) == 3:
        ### Code for connecting 3 parts 
        print("Only one input fdf for scattering region is required for now")
        sys.exit(111)

    ### generate fdf is not input
    if not lfdf:
        ### Use default side part in scattering region
        if not lsidepart:
            path_scatter_left = dir_default+'/elecL.fdf'
            path_scatter_rigt = dir_default+'/elecR.fdf'
        else:
            # Code for cli input
            pass

        elecL = modulename.read_struct(path_scatter_left)
        elecR = modulename.read_struct(path_scatter_rigt)

        
        if channel == 'grp':
            cnt = carbonlab.grp_rect(2,chsize)
            cnt.select_all()
            cnt.rotate(90, axis_dir=(1,0,0), with_cell=True)
            cnt.translate(-1.42, 0, 0)
            cnt.wrap_positions()
            cnt.select_z(-0.1, 1)
            cnt.delete()
            ### 2.2 Combine left and right parts of electrode
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
            fout = f"{channel}_{chsize}.fdf"
            path_fdf = cwd+f'/{fout}'
            modulename.write_struct(new, fname = f"{fout}")
            
            # write_poscar(new, f"cnt_{i}.poscar")
        else:
            ### different lab is required 
            pass

    ### 2.2 find psf files
    #print(f"in mdeling scattering region {new.get_species()}")
    atom_symbols = new.get_species()
    path_pp=[]
    for at in atom_symbols:
        f_psf = at+'.psf'
        path_psf = modulename.default_location+f'/psf/{f_psf}'
        path_pp.append(path_psf)
    
    '''
    if job == 'model':
        print(f"copy {path_fdf_left} to {cwd}")
        shutil.copy(path_fdf_left, cwd)
        shutil.copy(path_fdf_rigt, cwd)
        return 0
    '''
    
    dict_scatter={'psf': path_pp, 'fdf': path_fdf }

    return dict_elec, dict_scatter

def qtTransport(calc, dict_elec, dict_model, qt_dir, inp, outp, np=1, det_k=0, max_nk=10, init_nk=1, step_nk=1, opts=None):
    '''
    calc        instance of calculator such as Siesta
    dict_elec   dictionary of electrode {'psf': M.psf, 'fdf': [elec_left.fdf,elec_rigt.fdf]}
    dict_model  dictionary of scatter   {'psf': [M1.psf, M2.psf], 'fdf': scatter.fdf}
    qt_dir      [elec calc, scatter calc, postprocess]
    inp         yaml input file for scattering calculation
    outp        yaml output for scattering and input for postprocess
    '''
    ### 0. Models are generated in advance
    ### 1. Electrode calculation
    calcElectrode(calc, dict_elec, qt_dir[0], np)
    ### 2. Scattering region calculation
    calcScattering(calc, dict_model, qt_dir[1], inp, outp, np)
    ### 3. Post processing
    qtPlot(calc, qt_dir[2], qt_dir[1], outp)

    return 0

def calcElectrode(calc, el_model, elecdir, np):
    cwd = os.getcwd()
    calc.set_mode('elec')
    ### Make main subdirectory for electrode calculation
    if not os.path.exists(elecdir):
        os.mkdir(elecdir)
        os.chdir(elecdir)
        ### two calculations for left and right electrode
        ### make two subdirectories inside 1elec dir
        sub_dirs=[]
        f_psf = re.split('/', el_model['psf'])[-1]
        metal = os.path.splitext(f_psf)[0]
        sub_dirs.append(metal+'_left')
        sub_dirs.append(metal+'_rigt')
                
        for subdir, fstruct in zip(sub_dirs, el_model['fdf']):
            os.mkdir(subdir)
            os.chdir(subdir)
            os.mkdir('Input')
            shutil.copy(el_model['psf'], 'Input')
            shutil.copy(fstruct, 'Input/STRUCT.fdf')
            shutil.copytree('Input', 'Run')
            os.chdir('Run')
            calc.read_all_fdf()
            calc.runQuantumTransport(np)
            print(f"Job completed in {subdir}")
            os.chdir('..')
            os.chdir('..')
    os.chdir(cwd)
    return 0


def calcScattering(calc, dict_model, scatter_dir, finp, foutp, nproc):
    modulename = importlib.import_module(calc.__class__.__module__)
    '''
    finp    Additional input file of yaml
    '''
    dir_default = modulename.default_location
    cwd = os.getcwd()
    if not os.path.exists(scatter_dir):
        os.mkdir(scatter_dir)
    os.chdir(scatter_dir)
    ###### working inside work subdirectory
    cwdw = os.getcwd()
    if not os.path.isdir('Input'):
        os.mkdir('Input')
    ### copy model.fdf, elec.fdf, psfs
    shutil.copy(dict_model['fdf'], 'Input/STRUCT.fdf')
    for psf in dict_model['psf']:
        shutil.copy(psf, 'Input/')
        
    ### if not exist, copy yaml input file copy model from 2.model/
    if os.path.isfile(f"{cwd}/{finp}"):
        path_inp = f"{cwd}/{finp}"
    elif os.path.isfile(f"{dir_default}/{finp}"):
        path_inp = dir_default+'/'+finp
    else:
        print(f"Can't find {finp} file")
        sys.exit(121)
    shutil.copy(path_inp, '.')

    with open(finp) as f:
        yoption = yaml.safe_load(f)

    os.chdir('Input')
    calc.read_all_fdf()
    os.chdir('..')

    ### Develop for-loop in case voltage list, also mod format inside yaml
    voltage = yoption['Voltage'].split()[0]    # [1] is 'eV'
    #for voltage in voltages:
    if voltage not in os.listdir():
        vdir_name = voltage + 'eV'
        if not os.path.isdir(vdir_name):
            os.mkdir(vdir_name)
        os.chdir(vdir_name)
        cwdv = os.getcwd()          # voltage dir

        ### temparay modification not to run

        if not os.path.isdir('TSHS'):
            shutil.copytree('../Input', 'TSHS')
        os.chdir('TSHS')
        calc.set_mode('scatter')
        calc.set_option('TS.Elecs.Neglect.Principal', True)
        check_fname = f'{cwdv}/TSHS/MESSAGES'
        
        if (not os.path.isfile(check_fname)) or (not check_file(check_fname, job_complete)):
            calc.runQuantumTransport(nproc, **yoption)
        yoption['scatter'] = os.getcwd()
        os.chdir(cwdv)

        if not os.path.isdir('TBTrans'):
            shutil.copytree('../Input', 'TBTrans')
            os.chdir('TBTrans')
            calc.set_mode('tbtrans')
            calc.read_all_fdf()
            calc.set_option('kgrid_Monkhorst_Pack', True, (12,1,1))
            calc.set_option('TS.Elecs.Neglect.Principal', True)
            calc.runQuantumTransport(nproc, **yoption)
            yoption['tbtrans'] = os.getcwd()
        os.chdir(cwdw)
        ### write to yaml
        with open(f'{foutp}', 'w') as f:
            yaml.safe_dump(yoption, f)
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

