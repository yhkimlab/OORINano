###### Directory tree
### running files
run_qtnegf.py           main script
input.yaml              input parameters-dirname should be same with dirs in main script
Models/elec/            each Au electrode models are provided for semi-infinite electrode and 
      /channel/         - electrode parts in scattering region
### Auxiliary files
readme.txt, Makefile    for convinience
Ex_model/               directory to practice making graphene models
Test_params/            to check the change of parameters due to the user input parameters
                        filename should have 'el' or 'sc' for electrode and scattering calculation, respectively
cf. default parameters are $NCHOME/simulators/siesta/siesta_default/transmission/elec, scatter
    psf files locates on $NCHOME/simulators/siesta/siesta_default/psf

###### set env
###check $NCHOME/aux/env_slurm.py
    ln -s your_env.py env.py

###### Usage
options:
    -j [run (all calculation)|model (print model at workdir)|fdf (show fdf parameters)]            
    -c   atoms in scattering model            
    -cs   fdf scattering structure: 1 for channel, 2 for parts of electrode            
    -e   one atom for electrode            
    -es   fdf 2 electrode structure for left & right
E.G.::
    $ python run_qtnegf.py -j run   -c grp -cs 6 -e Au -jd 1.9 -np 20
    $ python run_qtnegf.py -j model -c grp -cs 6 -e Au -jd 1.9
    $ python run_qtnegf.py -j params -p Test_params/elec.fdf
    - before run python with -j run, make clean to remove subdirectories

###### Working procedure
1. make models for electrode and scattering region (model)
    model for scattering region is copied at work directory
2. electrode calculation in ./1elec
    input in    1elec/Au_left/Input
    siesta run  1elec/Au_left/Run
    default directory ~/siesta/siesta_default/transmission/elec
    if any fdf files in wdir, it will overwrite the key-values
3. transport calculation in ./2channel
    reference input file "input.yaml"
    TSHS and TBtrans runs in subdir of voltages such as 0.1eV/ 0.2eV/ 0.5eV/
    default directory ~/siesta/siesta_default/transmission/scatter
    if any fdf files in wdir, it will overwrite the key-values
    
    keyword in "output.yaml" is used for postprocessing
4. post_processing in ./3postprocess
    fdf2xcrysden.py result.jpg
    
    
### Practice generating models
/Example_models
 python generate_model.py




    

