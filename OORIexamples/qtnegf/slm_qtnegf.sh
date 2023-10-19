#!/bin/bash
#SBATCH --nodes=1
##
#SBATCH --time=90-12:34          # Runtime: Day-HH:MM
#SBATCH -o sms.%N.%j.o         # STDOUT 
#SBATCH -e sms.%N.%j.e         # STDERR
##

#module purge

#module add compiler/2022.1.0
#module add mkl/2022.1.0
#module add mpi/2021.6.0 

#export OMP_NUM_THREADS=1

pdir=$SLURM_SUBMIT_DIR
jobname=$SLURM_JOB_NAME
wdir=$pdir/$jobname
logfile=$pdir/${SLURM_JOB_ID}.${jobname}.log
outfile=$pdir/${SLURM_JOB_ID}.${jobname}.out
partname=$SLURM_JOB_PARTITION
nodelist=$SLURM_JOB_NODELIST

subdir=1    # change into other value to run main directory

str="-j run -c grp -cs 6 -e Au -jd 1.9 -np $SLURM_NTASKS -x $partname"

if [ $subdir -eq 1 ]; then

    if [ ! -f "$wdir" ]; then
        mkdir $wdir
    else
        echo "there exists $wdir"
    fi

    cp -r input.yaml Models $wdir
    cd $wdir
    
    exe="../run_qtnegf.py   $str"

else
    exe="run_qtnegf.py      $str"
fi

echo `date` >> $logfile
echo "Partition : $partname" >> $logfile
echo "python3 $exe" >> $logfile

python3 $exe >> $logfile

echo `date` >> $logfile

cd $pdir
mv $logfile $outfile
