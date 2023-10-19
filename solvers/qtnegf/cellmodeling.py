import sys, os, re, importlib
from ...models.interface import twoInterfaces

'''
modeling electrode
    if prepared, copy from simulator default directory
    if not, modeling from module models are required
'''

elec_atom = ['Au']

def get_psf(calc, files):
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
    else:
        print("Input error in ge_elec_inifile")
        sys.exit(112)

### separate path to model and psf
def model_electrode(calc, el_struct, path2model=None):
    '''
    calc    to obtain default directory pathway for psf
                      module readAtomicStructure
    To make electrode model
        el_struct       M or M.spf and/or 2 fdf file
    '''
    cwd = os.getcwd()
    
    simmodule = importlib.import_module(calc.__class__.__module__)
    
    ### extract fdf file
    li_fdf = [ f for f in el_struct if 'fdf' in f]
    ### remaining might be .psf
    res = list(filter(lambda x: x not in li_fdf, el_struct ))
    lfdf = False
    lpsf = False
    ### if Au_STRUCT_left[rigt].fdf exists
    if len(li_fdf) == 2:
        path_fdf_left = cwd+'/'+li_fdf[0]
        path_fdf_rigt = cwd+'/'+li_fdf[1]
        lfdf = True
        
    elif len(li_fdf) == 1:
        print("Can't run with one-side electrode structure")
        ### module for rotate left struct to right struct
        sys.exit(101)

    ### res might be M or M.psf
    if len(res) == 1:
        if os.path.splitext(res[0])[-1] == '.psf':
            path_psf = cwd+'/'+el_struct[0]
            metal = os.path.splitext(el_struct[0])[0]
            lpsf = True
        else:
            metal = el_struct[0]
            f_psf = metal + '.psf'
    ### No psf: read input file     
    elif len(res) == 0:
        metal = simmodule.readAtomicStructure(path_fdf_left)[0].get_symbol()
        f_psf = metal+'.psf'
    else:
        print(f"Algorithm error for electrode model input")
        sys.exit(103)
    
    ### if not path_psf, get default
    if not lpsf and 'f_psf' in locals():
        path_psf = get_psf(calc, f_psf)
    else:
        ### Codes for generation of electrode
        print("Electrode generation code is not ready")
        sys.exit(105)
    ### if no path_fdf, get default
    if not lfdf and 'metal' in locals():
        f_fdf1 = metal+'_STRUCT_left.fdf'
        f_fdf2 = metal+'_STRUCT_rigt.fdf'
        ### path to wdir/Models
        if path2model:
            #print(f"find model in {path2model}")
            path_fdf_left = f"{path2model}/{f_fdf1}"
            path_fdf_rigt = f"{path2model}/{f_fdf2}"
            if not os.path.exists(path_fdf_left) or not os.path.exists(path_fdf_rigt):
                print("files does not exist")
                sys.exit(103)
            ###  default dir in siesta module
            #path_fdf_left, path_fdf_rigt = get_elec_inifile(calc, [f_fdf1, f_fdf2])
            #print("define path to electrode models")

        ### Generation Codes are required 
        else:
            print("path to model is requisite")
            sys.exit(102)
   
    dict_elec={'psf': path_psf, 'fdf': [path_fdf_left, path_fdf_rigt]}
    return dict_elec

def model_channel(calc, channel_struct, chsize, junc_dist, path2model=None):
    '''
    calc    to obtain default directory pathway for psf
                      module readAtomicStructure
    To make scattering model
        3 or 1 fdf files
            2 side electrode parts
            1 scattering region
    '''
    cwd = os.getcwd()
    simmodule = importlib.import_module(calc.__class__.__module__)
    
    ###### get fdf file: [channel, left-elect, right-elect]
    ### if only channel given, get default electrode part
    #dir_default = simmodule.default_location+'/model/channel'
    
    lfdf = False
    lsidepart = False

    fdf = [f for f in channel_struct if 'fdf' in f]
    ### No inupt fdf but model name, gnenerate
    #print(f"{fdf}")
    ### total scattering region
    if len(fdf) == 1:
        path_fdf = f"{cwd}/{fdf[0]}"
        cell = simmodule.readAtomicStructure(path_fdf)
        lfdf = True
        
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

    ### Generate complete cell for channel
    if not lfdf:
        ### Use default side part in scattering region
        if not lsidepart:
            path_scatter_left = path2model+'/elecL.fdf'
            path_scatter_rigt = path2model+'/elecR.fdf'
            if not os.path.exists(path_scatter_left) or not os.path.exists(path_scatter_rigt):
                print(f"Can't find electrode parts for making channel model {path_scatter_rigt}")
                sys.exit(112)
        else:
            ### Codes for generation (electrode parts for channel) are required
            pass

        elecL = simmodule.readAtomicStructure(path_scatter_left)
        elecR = simmodule.readAtomicStructure(path_scatter_rigt)

        dict_channel = {'name': channel, 'size': chsize}
        cell = twoInterfaces(dict_channel, elecL, elecR, junc_dist)
        fout = f"{channel}_{chsize}.fdf"
        path_fdf = cwd+f'/{fout}'
        ### write channel in work directory
        simmodule.writeAtomicStructure(cell, fname = f"{fout}")
        lfdf = True
    
    ### find psf files: from input or default directory
    #print(f"in mdeling scattering region {new.get_species()}")
    path_pp=[]
    psf = [f for f in channel_struct if 'psf' in f]

    if len(psf) != 0:
        for pp in psf:
            path_psf = f"{cwd}/{pp}"
            path_pp.append(path_psf)
    else:
        atom_symbols = cell.get_species()
        
        for at in atom_symbols:
            f_psf = at+'.psf'
            path_psf = simmodule.default_location+f'/psf/{f_psf}'
            path_pp.append(path_psf)
    

    dict_scatter={'psf': path_pp, 'fdf': path_fdf }

    return dict_scatter

