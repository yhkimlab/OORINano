#!/bin/bash
#SBATCH -J testdos         # job name
#SBATCH -o stdout.txt      # output and error file name (%j expands to jobID)
#SBATCH -p X4 
#SBATCH -N 1               # total number of nodesmpi tasks requested
#SBATCH -n 24               # total number of mpi tasks requested

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
    exit 1
fi

runpython="run_catalysis.py"
poscar=${pos}
chgcar=${chg}
djob=spdos
rjob=${job:-$djob}
echo $rjob 
echo `date` > $logfile
### in case no CHGCAR, run sp and write CHGCAR
if [ $rjob == "spdos" ]; then
    echo $rjob is spdos
    cp $poscar $wdir/POSCAR
    cd $wdir
    str="../$runpython -j $rjob -p POSCAR -n $SLURM_JOB_NUM_NODES -np $SLURM_NTASKS --npar $npar"
    echo "python3 $str | vasp running" >> $logfile
    echo "python3  $str " >> $logfile
    python3 $str >> $logfile
    cd $pdir
    echo "sp was done" >> $logfile

    rjob=dos
    poscar=${wdir}/CONTCAR
    chgcar=${wdir}/CHGCAR
    wdir=${wdir}dos
    if [ ! -f "$wdir" ]; then
        mkdir $wdir
    else
        echo "there exists $wdir"
        exit 1
    fi
fi

cp $poscar $wdir/POSCAR
cp $chgcar $wdir
cd $wdir
str="../$runpython -j $rjob -p POSCAR -n $SLURM_JOB_NUM_NODES -np $SLURM_NTASKS --npar $npar"
echo "python3 $str | vasp running" >> $logfile
echo "python3  $str " >> $logfile
python3 $str >> $logfile
echo "dos calculation was done" >> $logfile
echo `date` >> $logfile
cd $pdir
mv $logfile $outfile

