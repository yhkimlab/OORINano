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

def qtNegf(calc, dict_elec, dict_model, qt_dir, inp, outp, np=1, det_k=0, max_nk=10, init_nk=1, step_nk=1, opts=None):
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
        print(f"{elecdir} was generated")
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
    simmodule = importlib.import_module(calc.__class__.__module__)
    '''
    finp    Additional input file of yaml
    '''
    dir_default = simmodule.default_location
    cwd = os.getcwd()
    if not os.path.exists(scatter_dir):
        os.mkdir(scatter_dir)
        print(f"{scatter_dir} was generated")
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

