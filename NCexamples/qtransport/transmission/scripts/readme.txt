###### setup env in $NCHOME/aux/env.py

###### Calculation of Quantum Transport
    direct run: python run_qt.py -n nproc

### Step 1: Electrode calculation in 1.elect
1. initial condition
    1.elect/Au_left[right]/input/  Au.psf, STRUCT.fdf
2. python ../1elec_calc.py -n nproc (run in 1.elect/)
    cp input to OUT
    write OUT/*

### Step 2: Scattering region modeling in 2.model
1. initial condition
    2.model/initial/  elecL.fdf, elecR.fdf
2. python ../2generate_model.py  (run in 2.model/)
    generate diverse sizes of zigzag-graphene models

### Step 3: Scattering region calculation in 3.scatter_tbtrans
1. initial condition
    input.yaml, input/psf files
    get one model in 2.model

2. python ../3scatter.py -n nproc  (run in 3.scatter_tbtrans)
    select one model 
        copy 2.model/anymodel to input/
    make subdirectory of voltage such as 0.1/ 0.2/ 0.5/
        TSHS/calculate scattering region with electrode and pseudo potential
            copy 1.elect/  electL.TSHS, elecR.TSHS
        TBtrans/calculate tight binding transport
    output: output.yaml
    
### Step 4: Post-processing in 4.post_processsing
1. initial condition: None
2. python ../4.post_process.py  (run in 4.post_processsing)

    
    




    

