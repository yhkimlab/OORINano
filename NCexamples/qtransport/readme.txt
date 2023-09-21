###### set env
###check $NCHOME/aux/env_slurm.py
    ln -s your_env.py env.py

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
    
    
### Test for model generation
/model
 python generate_model.py




    

