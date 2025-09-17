#!/bin/bash
#SBATCH -J testdos         # job name
#SBATCH -o stdout.txt      # output and error file name (%j expands to jobID)
#SBATCH -p X4 
#SBATCH -N 1               # total number of nodesmpi tasks requested
#SBATCH -n 24               # total number of mpi tasks requested

pdir=$SLURM_SUBMIT_DIR
jobname=$SLURM_JOB_NAME
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

runpython="run_vasp.py"
poscar="POSCAR."${jobname}

echo $vasjob >> $logfile
echo `date` >> $logfile

### Run sp at wdir :: in case no CHGCAR, run sp and write CHGCAR
wdir=${pdir}/${jobname}${vasjob}              # wdir is odir_name+'sp'
if [ $vasjob == "sp" -o $vasjob == "opt" -o $vasjob == "md"]; then
    ### make sp directory
    if [ ! -f "$wdir" ]; then
        mkdir $wdir
    else
        echo "there exists $wdir"
        exit 1
    fi
    cp $poscar $wdir/POSCAR
    cd $wdir
    str="../$runpython -j $vasjob -p POSCAR -n $SLURM_JOB_NUM_NODES -np $SLURM_NTASKS --npar $npar"
    echo "python3 $str | vasp running" >> $logfile
    echo "python3  $str " >> $logfile
    python3 $str >> $logfile
    cd $pdir
    echo "sp was done" >> $logfile

    poscar=${wdir}/CONTCAR
    chgcar=${wdir}/CHGCAR
else
    wdir=${pdir}/${old_dir}     # old_dir from local
fi


if [ $vasjob == "md" ]; then
    

echo `date` >> $logfile
cd $pdir
mv $logfile $outfile
