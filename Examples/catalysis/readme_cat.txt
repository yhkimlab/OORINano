### Order for test
1. download NanoCore and change dirname to nanocore
    git clone https://github.com/yhkimlab/NanoCore nanocore
2. set PYTHONPATH upto parent directory of nanocore
3. copy ECexamples to your working directory
4. run job script as follows

### Details for Running

1. Set VASP env
OORIHOME = OORINano home directory
$OORIHOME/utils/env.py or your_env.py
    e.g.
        $ln -s your_env.py env.py
    set vasp_calculator -> to your vasp binary file
    set vasp_POTCAR_dft -> to your POTCAR directory

2. How to run
'orr'   ORR (Oxygen   Reduction Reaction)
'her'   HER (Hydrogen Evolution Reaction)
also check run_catalysis.py
    $python run_catalysis.py -u

a. slurm with partition
    sbatch -J testorr -p X5 -n 1 -n 32 [--export=cat='orr'] [--export=pos='cp'] slurm_sbatch_nc.sh
    sbatch -J testher -p X5 -n 1 -n 32  --export=cat='her'   --export=pos='gen' slurm_sbatch_nc.sh
        qname       -J      make dir in the qname and vasp runs in the dir       
        partition   -p
        nNode       -n
        nproc       -np
        [default]   cat='orr' ['her']
                    pos='cp'  anyword: 'cp' copies prepared POSCAR, otherwise, provide -p to generate metal slab
b. pbs without partition
    qsub -N testorr [-v cat='orr'] pbs_vasp_kisti_skl.sh
c. direct run
    run_catalysis.py -j orr -sj run -n 1 -np 24 [--npar $npar|--ncore $ncore] 
