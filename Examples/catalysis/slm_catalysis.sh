#!/bin/bash
#SBATCH -J TEST-NC         # job name
#SBATCH -o stdout.txt      # output and error file name (%j expands to jobID)
#SBATCH -p X4 
#SBATCH -N 1               # total number of nodesmpi tasks requested
#SBATCH -n 24               # total number of mpi tasks requested

###### Usage
### SLURM: $sbatch -J test -p X3 -N 4 -n 80 slurm_sbatch_nc.sh
### Job runs in "test" dir
### .log turns into .out when job finished

pdir=$SLURM_SUBMIT_DIR
jobname=$SLURM_JOB_NAME
wdir=$pdir/$jobname
logfile=$pdir/${SLURM_JOB_ID}.${jobname}.log
outfile=$pdir/${SLURM_JOB_ID}.${jobname}.out
partname=$SLURM_JOB_PARTITION
nodelist=$SLURM_JOB_NODELIST

if [ ${partname} == 'X1' ] || [ ${partname} == 'X1c' ]; then
    par=2; SLURM_CPUS_PER_NODE=8
elif [ ${partname} == 'X2' ] || [ ${partname} == 'X2c']; then
    par=2; SLURM_CPUS_PER_NODE=12
elif [ $partname == 'X3' ]; then
    par=4; SLURM_CPUS_PER_NODE=20
elif [ $partname == 'X4' ]; then
    par=4; SLURM_CPUS_PER_NODE=24
else    # if X5 or X6
    par=4; SLURM_CPUS_PER_NODE=32
fi

npar=$(expr $SLURM_JOB_NUM_NODES \* $par )


if [ ! -d "$wdir" ]; then
    mkdir $wdir
else
    echo "there exists $wdir"
fi

jobfile="run_catalysis.py"
catkind=${cat:-"orr"}     # select [ORR|HER]
poscar=${pos:-"gen"}

if [ $poscar == 'gen' ]; then
    if [ $catkind == 'orr' ]; then
        str="../$jobfile -j run -r orr -m Pt -ss 111 3 -n $SLURM_JOB_NUM_NODES -np $SLURM_NTASKS --npar $npar "
    else # HER
        str="../$jobfile -j run -r $catkind -m Pt Au Ag Pd Ni Cu -ss 111 3 -n $SLURM_JOB_NUM_NODES -np $SLURM_NTASKS --npar $npar "
    fi
else
    cp CONTCAR_Pt-SAC $wdir/POSCAR
    str="../$jobfile -j run -r $catkind -i POSCAR -n $SLURM_JOB_NUM_NODES -np $SLURM_NTASKS --npar $npar"
fi

echo `date` >> $logfile
echo HOME is $HOME >> $logfile
cd $wdir
echo " Partition : $SLURM_JOB_PARTITION" >> $logfile
echo " Nodelist  : $nodelist " >> $logfile
echo "python3 $jobfile running | vasp running" >> $logfile
echo "python3  $str " >> $logfile
python3 $str >> $logfile
echo `date` >> $logfile

cd $pdir
mv $logfile $outfile
