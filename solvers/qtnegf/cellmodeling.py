import sys, os, re, importlib
from ...models.interface import twoInterfaces

'''
modeling electrode
    if prepared, copy from simulator default directory
    if not, modeling from module models are required
'''

elec_atom = ['Au']

def get_elec_inifile(calc, files):
    '''
    to get default pathway, calc is passed
    '''
    cwd = os.getcwd()
    simmodule = importlib.import_module(calc.__class__.__module__)
    
    ### .psf file
    if type(files) == str and re.search('psf', files):
        dir_default = simmodule.default_location+'/psf'
        fpath = dir_default+'/'+files
        if os.path.isfile(fpath):
            return fpath
        else:
            print(f"Can't find {fpath}")
            sys.exit(111)
    ### electrode 2 fdf files
    elif type(files) == list:
        dir_default = simmodule.default_location+'/model/elec'
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

def model_electrode(calc, el_struct, el_size):
    '''
    calc                to obtain default directory pathway
    To make electrode model
        el_struct       M or M.spf and/or 2 fdf file
    '''
    cwd = os.getcwd()
    
    simmodule = importlib.import_module(calc.__class__.__module__)
    
    li_fdf = [ f for f in el_struct if 'fdf' in f]
    lfdf = False
    lpsf = False
    ### if Au_STRUCT_left[rigt].fdf exists
    if len(li_fdf) == 2:
        path_fdf_left = cwd+'/'+li_fdf[0]
        path_fdf_rigt = cwd+'/'+li_fdf[1]
        lfdf = True
        el_struct = [ f for f in el_struct if f not in fdf ]
    elif len(li_fdf) == 1:
        print("Can't run with one-side electrode structure")
        ### module for rotate left struct to right struct
        sys.exit(101)

    ### if metal species exists
    #li_psf = [ f for f in el_struct if 'psf' in f]
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
        metal = simmodule.readAtomicStructure(path_fdf_left)[0].get_symbol()
        f_psf = metal+'.psf'
    
    ### if not path_psf, get default
    if not lpsf and 'f_psf' in locals():
        path_psf = get_elec_inifile(calc, f_psf)
    else:
        ### Codes for generation of electrode
        print("Electrode generation code is not ready")
        sys.exit(102)
    ### if no path_fdf, get default
    if not lfdf and 'metal' in locals():
        f_fdf1 = metal+'_STRUCT_left.fdf'
        f_fdf2 = metal+'_STRUCT_rigt.fdf'
        path_fdf_left, path_fdf_rigt = get_elec_inifile(calc, [f_fdf1, f_fdf2])
    
    dict_elec={'psf': path_psf, 'fdf': [path_fdf_left, path_fdf_rigt]}
    return dict_elec

def model_channel(calc, channel_struct, chsize, junc_dist):
    cwd = os.getcwd()
    simmodule = importlib.import_module(calc.__class__.__module__)
    
    ###### get fdf file: [channel, left-elect, right-elect]
    ### if only channel given, get default electrode part
    dir_default = simmodule.default_location+'/model/channel'
    
    lfdf = False
    lsidepart = False

    fdf = [f for f in channel_struct if 'fdf' in f]
    ### No inupt fdf but model name, gnenerate
    #print(f"{fdf}")
    if len(fdf) == 1:
        lfdf = True
        pass
    elif len(fdf) == 0:
        ### if input is one, it is channel name
        if len(channel_struct) == 1:
            channel = channel_struct[0]
        else:
            print(f"Algorithm error in preparing scattering region")
            sys.exit(110)
    elif len(fdf) == 3:
        ### Code for connecting 3 parts 
        print("Only one input fdf for scattering region is required for now")
        sys.exit(111)

    ### Generate channel model
    if not lfdf:
        ### Use default side part in scattering region
        if not lsidepart:
            path_scatter_left = dir_default+'/elecL.fdf'
            path_scatter_rigt = dir_default+'/elecR.fdf'
        else:
            # Code for cli input
            pass

        elecL = simmodule.readAtomicStructure(path_scatter_left)
        elecR = simmodule.readAtomicStructure(path_scatter_rigt)

        dict_channel = {'name': channel, 'size': chsize}
        cell = twoInterfaces(dict_channel, elecL, elecR, junc_dist)
        fout = f"{channel}_{chsize}.fdf"
        path_fdf = cwd+f'/{fout}'
        simmodule.writeAtomicStructure(cell, fname = f"{fout}")
        lfdf = True
    
    ### find psf files
    #print(f"in mdeling scattering region {new.get_species()}")
    atom_symbols = cell.get_species()
    path_pp=[]
    for at in atom_symbols:
        f_psf = at+'.psf'
        path_psf = simmodule.default_location+f'/psf/{f_psf}'
        path_pp.append(path_psf)
    
    '''
    if job == 'model':
        print(f"copy {path_fdf_left} to {cwd}")
        shutil.copy(path_fdf_left, cwd)
        shutil.copy(path_fdf_rigt, cwd)
        return 0
    '''
    
    dict_scatter={'psf': path_pp, 'fdf': path_fdf }

    return dict_scatter

