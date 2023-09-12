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

job=HER     # select [ORR|HER]

pdir=$SLURM_SUBMIT_DIR
jobname=$SLURM_JOB_NAME
wdir=$pdir/$jobname
logfile=$pdir/${SLURM_JOB_ID}.${jobname}.log
outfile=$pdir/${SLURM_JOB_ID}.${jobname}.out

if [ ! -f "$wdir" ]; then
    mkdir $wdir
else
    echo "there exists $wdir"
    exit 1
fi

### copy files to $wdir
if [ $job == "ORR" ]; then
    cp test_ORR.py $wdir
    cp CONTCAR_Pt-SAC $wdir/POSCAR
elif [ $job == "HER" ]; then
    cp test_HER.py $wdir
fi

echo `date` > $logfile
cd $wdir

python3 ./test_HER.py
echo `date` >> $logfile

cd $pdir
mv $logfile $outfile
