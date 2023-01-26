###### How to run ORR, HER using sbatch
### ORR
sbatch -J orr -p X5 -N 1 -n 32 --export=job='orr' slurm_sbatch_nc.sh
### HER
sbatch -J her -p X5 -N 1 -n 32 --export=job='her' slurm_sbatch_nc.sh
### qname: -J
### partition: -p
### nNode: -N
### nproc: -n
