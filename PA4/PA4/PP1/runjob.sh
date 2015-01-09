#
# Batch submission script for Prog. Assignment 2
# If this script is placed in file runjob.sh,
# to submit batch job, use "qsub runjob.sh"
# Status of batch job can be checked using "qstat -u `whoami`" 
#
#PBS -l walltime=0:30:00
#PBS -l nodes=64:ppn=1
#PBS -N PA4
#PBS -S /bin/ksh
#PBS -j oe

INP=$HOME/PA4
#TMPDIR=`mktemp -d`
pbsdcp $INP/a.out $TMPDIR
cd $TMPDIR
echo " "; echo "--------------------------------------------------------------"
echo "Programming Assignment 4, Problem 1"
echo "--------------------------------------------------------------"; echo " "
#mpicc sieve.c 
#pbsdcp ./a.out $TMPDIR
#cd $TMPDIR
echo " "; echo "--- 1 core ---"
mpiexec -n 1 ./a.out 1000000000
echo " "; echo "--- 2 core ---"
mpiexec -n 2 ./a.out 1000000000
echo " "; echo "--- 4 core ---"
mpiexec -n 4 ./a.out 1000000000
echo " "; echo "--- 8 core ---"
mpiexec -n 8 ./a.out 1000000000
echo " "; echo "--- 16 core ---"
mpiexec -n 16 ./a.out 1000000000
echo " "; echo "--- 32 core ---"
mpiexec -n 32 ./a.out 1000000000
echo " "; echo "--- 64 core ---"
mpiexec -n 64 ./a.out 1000000000
#echo " "; echo "--------------------------------------------------------------"
#echo "Programming Assignment 3, Problem 2"
#echo "--------------------------------------------------------------"; echo " "
#nvcc -O3 -arch=sm_20 PA3-Prob2.cu; ./a.out

rm -rf $TMPDIR
