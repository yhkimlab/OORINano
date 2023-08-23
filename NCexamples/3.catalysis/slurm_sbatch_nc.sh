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

if [ $partname == 'X1' ]; then
    par=2; SLURM_CPUS_PER_NODE=8
elif [ $partname == 'X2' ]; then
    par=2; SLURM_CPUS_PER_NODE=12
elif [ $partname == 'X3' ]; then
    par=4; SLURM_CPUS_PER_NODE=20
    #if [ $hmem -eq 1 ]; then
    if [ $hmem ]; then
        par=2
    fi
elif [ $partname == 'X4' ]; then
    par=4; SLURM_CPUS_PER_NODE=24
    if [ $hmem  ]; then
        par=2
    fi
else    # if X5
    par=4; SLURM_CPUS_PER_NODE=32
fi

npar=$(expr $SLURM_JOB_NUM_NODES \* $par )

if [ ! -f "$wdir" ]; then
    mkdir $wdir
else
    echo "there exists $wdir"
#    exit 1
fi

cat='orr'
cat_kind=${catkind:-$cat}     # select [ORR|HER]

jobfile="run_catalysis.py"
### copy files to $wdir
#if [ $cat_kind == 'orr' ]; then
#    cp CONTCAR_Pt-SAC $wdir/POSCAR
#fi

echo `date` > $logfile
cd $wdir
echo " Partition : $SLURM_JOB_PARTITION" >> $logfile
echo "python3 $jobfile running | vasp running" >> $logfile
#str="../$jobfile -j run -c $cat_kind -N $SLURM_JOB_NUM_NODES -np $SLURM_NTASKS --npar $npar"
str="../$jobfile -j run -c $cat_kind -p Pt 111 3 -N $SLURM_JOB_NUM_NODES -np $SLURM_NTASKS --npar $npar"
echo "python3  $str " >> $logfile
python3 $str >> $logfile
echo `date` >> $logfile

cd $pdir
mv $logfile $outfile
