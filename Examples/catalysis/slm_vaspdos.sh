#!/bin/bash
#SBATCH -J testdos         # job name
#SBATCH -o stdout.txt      # output and error file name (%j expands to jobID)
#SBATCH -p X4 
#SBATCH -N 1               # total number of nodesmpi tasks requested
#SBATCH -n 24               # total number of mpi tasks requested

pdir=$SLURM_SUBMIT_DIR
jobname=$SLURM_JOB_NAME
odir=$pdir/$jobname         # it is old directory of catalysis cal
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

runpython="run_vaspdos.py"
poscar=${pos}
chgcar=${chg}
djob=dos
rjob=${job:-$djob}
echo $rjob >> $logfile
echo `date` >> $logfile

### in case no CHGCAR, run sp and write CHGCAR
if [[ $rjob == *"sp"* ]]; then
    ### make sp directory
    wdir=${odir}sp              # wdir is odir_name+'sp'
    if [ ! -f "$wdir" ]; then
        mkdir $wdir
    else
        echo "there exists $wdir"
        exit 1
    fi
    echo $rjob has 'sp'
    cp $poscar $wdir/POSCAR
    cd $wdir
    str="../$runpython -j sp -p POSCAR -n $SLURM_JOB_NUM_NODES -np $SLURM_NTASKS --npar $npar"
    echo "python3 $str | vasp running" >> $logfile
    echo "python3  $str " >> $logfile
    python3 $str >> $logfile
    cd $pdir
    echo "sp was done" >> $logfile

    rjob=dos
    poscar=${wdir}/CONTCAR
    chgcar=${wdir}/CHGCAR
fi

if [[ $rjob == *"dos"* ]]; then
    wdir=${wdir}dos             # wdir is odir_name [+ 'sp'] + 'dos'
    if [ ! -f "$wdir" ]; then
        mkdir $wdir
    else
        echo "there exists $wdir"
        exit 1
    fi

    cp $poscar $wdir/POSCAR
    cp $chgcar $wdir
    cd $wdir
    str="../$runpython -j dos -p POSCAR -n $SLURM_JOB_NUM_NODES -np $SLURM_NTASKS --npar $npar"
    echo "python3 $str | vasp running" >> $logfile
    echo "python3  $str " >> $logfile
    python3 $str >> $logfile
    echo "dos calculation was done" >> $logfile
fi

echo `date` >> $logfile
cd $pdir
mv $logfile $outfile
