### 1.elect
1. initial condition
    input/Au.psf, STRUCT.fdf
2. (run) $ python 1.elec_calc.py -n nproc (direct run)
    cp input to OUT
    write OUT/*

### 2.model
1. initial condition
    initial/elecL.fdf, elecR.fdf
2. (run) $ python 2.generate_model.py
    generate diverse sizes of zigzag-graphene models

### 3. scatter_tbtrans
1. initial condition
    input.yaml, input/psf files
    one model in 2.model

2. (run) python 3.scatter.py -n nproc
    select one model 
        copy 2.model/anymodel to input/
    make subdirectory of voltage such as 0.1/ 0.2/ 0.5/
        TSHS/calculate whole region with pseudo potential
            electL.TSHS, elecR.TSHS, etc
        TBtrans/calculate tight binding transport
    output: output.yaml
    
### 4. post process
1. initial condition: None
2. (run) python 4.post_process.py
    
    




    

