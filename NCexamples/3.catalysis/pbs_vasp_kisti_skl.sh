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
ncpackage=$HOME/NanoCore/nanocore/
env_dir=$ncpackage/etc
rm $env_dir/env.py
ln -s $env_dir/env_kisti.py env.py

echo $jobname > $logfile
NPROC=`wc -l < $PBS_NODEFILE`
echo "NPROC = $NPROC" >> $logfile
echo start >> $logfile
date >> $logfile

cat='orr'
cat_kind=${catkind:-$cat}     # select [ORR|HER]

jobfile="run_catalysis.py"
### copy files to $wdir
#if [ $cat_kind == 'orr' ]; then
#    cp CONTCAR_Pt-SAC $wdir/POSCAR
#fi

if [ ! -f "$wdir" ]; then
    mkdir $wdir
else
    echo "there exists $wdir"
#    exit 1
fi


cd $log_dir/$wdir

str="../$jobfile -j run -c $cat_kind -p Pt 111 3 -np $NPROC --ncore 10"
echo "python3  $str " >> $logfile
python3 $str >> $logfile
echo `date` >> $logfile
echo "end" >> $logfile

mv $logfile $jobfile


