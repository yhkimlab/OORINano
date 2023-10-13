from __future__ import print_function
import sys, os, re, shutil, importlib
from ...aux import check_file

import yaml             ### added for cal scattering
from .qtplot import qtPlot

"""
Module QT-NEGF

Module for representing a quantum transport calculation

Parameters
----------
    dict_elect      {'psf': path2psf, 'fdf': path2fdf}
    dict_channel    {'psf': path2psf, 'fdf': path2fdf}
    qt_dir          dirnames for calc of electrode, channel, & post-process
"""

job_complete = "Job completed"

def qtNegf(calc, dict_elec, dict_model, qt_dir, inp, outp, fdf_params, np=1, det_k=0, max_nk=10, init_nk=1, step_nk=1, opts=None, show_params=False):
    '''
    calc        instance of calculator such as Siesta
    dict_elec   dictionary of electrode {'psf': M.psf, 'fdf': [elec_left.fdf,elec_rigt.fdf]}
    dict_model  dictionary of scatter   {'psf': [M1.psf, M2.psf], 'fdf': scatter.fdf}
    qt_dir      [elec calc, scatter calc, postprocess]
    inp         yaml input file for scattering calculation
    outp        yaml output for scattering and input for postprocess
    fdf_params  input parameters from work directory
                ['elec.fdf' (for electrode), 'scatter.fdf' (for scattering)]
    '''
    ### 0. Models are generated in advance
    ### 1. Electrode calculation
    fdf_elec    = None
    fdf_scatt   = None
    if fdf_params:
        for f in fdf_params:
            if re.search('el', f):
                fdf_elec = f
            elif re.search('sc', f):
                fdf_scatt = f
    calcElectrode(calc, dict_elec, qt_dir[0], np, fdf_elec = fdf_elec, show_params = show_params)
    ### 2. Scattering region calculation
    ### need to initialize Siesta instance to make data default
    calc = calc.__class__()
    calcScattering(calc, dict_model, qt_dir[1], inp, outp, np, fdf_scatt = fdf_scatt, show_params = show_params)
    if show_params:
        return 0
    else:
        ### 3. Post processing
        qtPlot(calc, qt_dir[2], qt_dir[1], outp)
        return 0

def calcElectrode(calc, el_model, elecdir, np, fdf_elec = None, show_params=False):
    cwd = os.getcwd()
    calc.set_mode('elec')       # sets default fdf params
    if fdf_elec:
        calc.read_fdf(fdf_elec)
    if show_params:
        print("=== Params for Electrode Calculation ===")
        calc.print_fdfs()
        return 0 
    ### Make main subdirectory for electrode calculation
    if not os.path.exists(elecdir):
        os.mkdir(elecdir)
        print(f"{elecdir} was generated")
    os.chdir(elecdir)
    ### two calculations for left and right electrode
    ### make two subdirectories inside elec subdir
    sub_dirs=[]
    f_psf = re.split('/', el_model['psf'])[-1]
    metal = os.path.splitext(f_psf)[0]
    sub_dirs.append(metal+'_left')
    sub_dirs.append(metal+'_rigt')
            
    for subdir, fstruct in zip(sub_dirs, el_model['fdf']):
        if not os.path.exists(subdir):
            os.mkdir(subdir)
        os.chdir(subdir)
        if not os.path.exists('Input'):
            os.mkdir('Input')
        shutil.copy(el_model['psf'], 'Input')
        shutil.copy(fstruct, 'Input/STRUCT.fdf')
        shutil.copytree('Input', 'Run')
        os.chdir('Run')
        cwdr = os.getcwd()
        check_fname = f'{cwdr}/MESSAGES'
        calc.read_all_fdf() # read all .fdf in /RUN
        if (not os.path.isfile(check_fname)) or (not check_file(check_fname, job_complete)):
            calc.runQtNegf(np)
            print(f"Job completed in {subdir}")
        os.chdir('..')
        os.chdir('..')
    os.chdir(cwd)
    return 0


def calcScattering(calc, dict_model, scatter_dir, finp, foutp, np, fdf_scatt = None, show_params=False):
    simmodule = importlib.import_module(calc.__class__.__module__)
    '''
    finp    Additional input file of yaml
    '''
    dir_default = simmodule.default_location
    cwd = os.getcwd()
    if not os.path.exists(scatter_dir):
        os.mkdir(scatter_dir)
    os.chdir(scatter_dir)
    ###### cwds: Scattering subdirectory
    cwds = os.getcwd()
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
    calc.read_all_fdf() # read struct
    os.chdir('..')

    ### Develop for-loop in case voltage list, also mod format inside yaml
    voltage = yoption['Voltage'].split()[0]    # [1] is 'eV'
    #for voltage in voltages:
    if voltage not in os.listdir():
        vdir_name = voltage + 'eV'
        if not os.path.isdir(vdir_name):
            os.mkdir(vdir_name)
        os.chdir(vdir_name)
        ### Voltage subdirectory
        cwdv = os.getcwd()

        ### temparay modification not to run

        if not os.path.isdir('TSHS'):
            shutil.copytree('../Input', 'TSHS')
        os.chdir('TSHS')
        calc.set_mode('scatter')
        ### Overwrite input params
        if fdf_scatt:
            calc.read_fdf(f"{cwd}/{fdf_scatt}")
        if show_params:
            print("=== Params for Scattering Calculation ===")
            calc.print_fdfs()
            return 0
        calc.set_option('TS.Elecs.Neglect.Principal', True)
        check_fname = f'{cwdv}/TSHS/MESSAGES'
        
        if (not os.path.isfile(check_fname)) or (not check_file(check_fname, job_complete)):
            calc.runQtNegf(np, **yoption)
        yoption['scatter'] = os.getcwd()
        os.chdir(cwdv)

        if not os.path.isdir('TBTrans'):
            shutil.copytree('../Input', 'TBTrans')
            os.chdir('TBTrans')
            calc.set_mode('tbtrans')
            calc.read_all_fdf()
            calc.set_option('kgrid_Monkhorst_Pack', True, (12,1,1))
            calc.set_option('TS.Elecs.Neglect.Principal', True)
            calc.runQtNegf(np, **yoption)
            yoption['tbtrans'] = os.getcwd()
        os.chdir(cwds)
        ### write to yaml
        with open(f'{foutp}', 'w') as f:
            yaml.safe_dump(yoption, f)
    os.chdir(cwd)
        
    return 0

