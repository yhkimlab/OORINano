#!/bin/sh
#PBS -V
#PBS -A vasp
#PBS -q normal
#PBS -l select=10:ncpus=40:mpiprocs=40:ompthreads=1
#PBS -l walltime=48:00:00

if [ -z $PBS_JOBNAME ]; then
    echo "Usage:: qsub -N dirname $SB/pypbs/pbs_vasp.sh"
    exit 1
fi

workdir=$PBS_O_WORKDIR
jobname=$PBS_JOBNAME
jobdir=$workdir/$jobname
logfile=$workdir/$jobname.log
outfile=$workdir/$jobname.out

echo $jobname > $logfile
echo $jobdir >> $logfile
NPROC=`wc -l < $PBS_NODEFILE`
echo "NPROC = $NPROC" >> $logfile
echo start >> $logfile
date >> $logfile

if [ ! -f "$jobdir" ]; then
    mkdir $jobdir
else
    echo "there exists $jobdir"
fi

jobfile="run_catalysis.py"
catkind=${cat:-"orr"}     # select [ORR|HER]
poscar=${pos:-"gen"}

if [ $poscar == 'gen' ]; then
	str="../$jobfile -j run -r $catkind -m Pt 111 3 -np $NPROC --ncore 10"
else
    cp $workdir/$pos $jobdir/POSCAR
	str="../$jobfile -j run -r $catkind -i POSCAR -np $NPROC --ncore 10"
fi

cd $jobdir

echo "python3  $str " >> $logfile
python3 $str >> $logfile
echo `date` >> $logfile
echo "end" >> $logfile

mv $logfile $outfile
