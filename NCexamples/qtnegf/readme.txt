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
    $ python run_qtnegf.py -j run -c grp -cs 5 -e Au -jd 1.9 -np 20
    $ python run_qtnegf.py -j model -c grp -cs 5 -e Au -jd 1.9

###### Working procedure
1. make models for electrode and scattering region (model)
    model for scattering region is copied at work directory
2. electrode calculation in ./1elec
    input in    1elec/Au_left/Input
    siesta run  1elec/Au_left/Run
3. transport calculation in ./2channel
    reference input file "input.yaml"
    TSHS and TBtrans runs in subdir of voltages such as 0.1eV/ 0.2eV/ 0.5eV/
    
    keyword in "output.yaml" is used for postprocessing
4. post_processing in ./3postprocess
    fdf2xcrysden.py result.jpg
    
    
### Practice generating models
/model
 python generate_model.py




    

