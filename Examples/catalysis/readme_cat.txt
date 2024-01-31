### Order for test
1. download OORINano and change dirname to oorinano
    git clone https://github.com/yhkimlab/OORINano oorinano
2. set PYTHONPATH upto parent directory of oorinano
3. copy ECexamples to your working directory
4. run job script as follows

### Details for Running

1. Set VASP env
OORIHOME = oorinano home directory
$OORIHOME/utils/env.py or your_env.py
    e.g.
        $ln -s your_env.py env.py
    set vasp_calculator -> to your vasp binary file
    set vasp_POTCAR_dft -> to your POTCAR directory

2. How to run catalysis
'orr'   ORR (Oxygen   Reduction Reaction)
'her'   HER (Hydrogen Evolution Reaction)
'oer'   OER (Oxygen   Evolution Reaction)
also check run_catalysis.py
    $python run_catalysis.py -u

a. slurm with partition
    sbatch -J testorr -p X5 -N 1 -n 32 [--export=cat='orr'] [--export=pos='cp'] slm_catalysis.sh
    sbatch -J testher -p X5 -N 1 -n 32  --export=cat='her',pos='gen' slm_catalysis.sh
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

3. Postprocess: DOS calculation
    After catalysis calculation
    3.1 sbatch

    3.2 Direct run
        3.2a Make a new dir and copy a poscar and run sp to make CHGCAR
            $python ../run_catalysis.py -j spdos -np 24 --npar 4
        3.2b In a new dir, run dos calculation with large k-points
            $python ../run_catalysis.py -j dos -np 24 --npar 4
