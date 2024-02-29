#!/bin/sh
#PBS -V
#PBS -A vasp
#PBS -q normal
#PBS -l select=10:ncpus=40:mpiprocs=40:ompthreads=1
#PBS -l walltime=24:00:00

if [ -z $PBS_JOBNAME ]; then
    echo "Usage:: qsub -N dirname $SB/pypbs/pbs_vasp.sh"
    exit 1
fi

log_dir=$PBS_O_WORKDIR
jobname=$PBS_JOBNAME
wdir=$jobname
logfile=$log_dir/$jobname.log
outfile=$log_dir/$jobname.out

### set env file here
ncpackage=$HOME/NanoCore/oorinano/
env_dir=$ncpackage/etc
rm $env_dir/env.py
ln -s $env_dir/env_kisti.py env.py

echo $jobname > $logfile
NPROC=`wc -l < $PBS_NODEFILE`
echo "NPROC = $NPROC" >> $logfile
echo start >> $logfile
date >> $logfile

if [ ! -f "$wdir" ]; then
    mkdir $wdir
else
    echo "there exists $wdir"
fi

jobfile="run_catalysis.py"
catkind=${cat:-"orr"}     # select [ORR|HER]
poscar=${pos:-"cp"}

if [ $poscar == 'cp' ]; then
    cp CONTCAR_Pt-SAC $wdir/POSCAR
	str="../$jobfile -j run -c $catkind -p POSCAR -np $NPROC --ncore 10"
else
	str="../$jobfile -j run -c $catkind -p Pt 111 3 -np $NPROC --ncore 10"
fi

cd $log_dir/$wdir

echo "python3  $str " >> $logfile
python3 $str >> $logfile
echo `date` >> $logfile
echo "end" >> $logfile

mv $logfile $jobfile


