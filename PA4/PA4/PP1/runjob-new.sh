#
# Batch submission script for Prog. Assignment 4
# If this script is placed in file runjob.sh,
# to submit batch job, use "qsub runjob.sh"
# Status of batch job can be checked using "qstat -u `whoami`" 
#
# TODO: Change nodes to be the maximum number of nodes you want to allocate
#
#PBS -l walltime=0:06:00
#PBS -l nodes=6:ppn=12
#PBS -N PA4-new
#PBS -S /bin/ksh
#PBS -j oe

TMPDIR=`mktemp -d --tmpdir=.`
cd $TMPDIR

#
# Assumes the programs are in a
# subdirectory called pa4 in your OSC home directory
# Change to match your directory and file name
#
cp $HOME/PA4/sieve.c .
mpicc -O3 sieve.c
#cp $HOME/pa4/MyMPI.* .
#mpicc -O3 sieve.c MyMPI.c 

for PROCS in 1 2 4 8 16 32 64;
do
    echo " "
    echo "--------------------------------------------------------------"
    echo "Programming Assignment 4, Sieve, $PROCS processes"
    echo "--------------------------------------------------------------"
    echo " "

    mpiexec -n $PROCS ./a.out 1000000000
done

rm -rf $TMPDIR
