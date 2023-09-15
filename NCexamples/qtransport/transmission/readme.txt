###### set env
    check $NCHOME/aux/env.py
    ln -s your_env.py env.py

###### Followings are run in one command consecutively
    refer to run_qt.py: $ python run_qt.py -u
### 1elec_cal.py run in 1.elec/
1. initial condition
    input/Au.psf, STRUCT.fdf
2. (run) $ python 1.elec_calc.py -np nproc (direct run)
    cp input to OUT
    write OUT/*

### 2generate_model.py runs in 2.model/
1. initial condition
    initial/elecL.fdf, elecR.fdf
2. (run) $ python 2generate_model.py
    generate diverse sizes of zigzag-graphene models

### 3scatter.py runs in 3.scatter_tbtrans/
1. initial condition
    input.yaml, input/psf files
    one model in 2.model

2. (run) python 3scatter.py -np nproc
    select one model 
        copy 2.model/anymodel to input/
    make subdirectory of voltage such as 0.1/ 0.2/ 0.5/
        TSHS/calculate whole region with pseudo potential
            electL.TSHS, elecR.TSHS, etc
        TBtrans/calculate tight binding transport
    output: output.yaml
    
### 4post_process.py runs in 4.post_processing/
1. initial condition: None
2. (run) python 4post_process.py
    
    




    

