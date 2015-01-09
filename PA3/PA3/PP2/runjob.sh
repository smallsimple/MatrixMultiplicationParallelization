#
# Batch submission script for Prog. Assignment 2
# If this script is placed in file runjob.sh,
# to submit batch job, use "qsub runjob.sh"
# Status of batch job can be checked using "qstat -u `whoami`" 
#
#PBS -l walltime=0:10:00
#PBS -l nodes=1:ppn=12:gpus=1
#PBS -N PA3
#PBS -S /bin/ksh
#PBS -j oe

module load cuda

TMPDIR=`mktemp -d`
cd $TMPDIR

#
# Assumes the cuda programs are in files PA3-Prob{1,2}.cu in a
# subdirectory called pa3 in your OSC home directory
# Change to match your directory and file name
#
cp $HOME/PA3/PP2/PA3-Prob2.cu .
#cp $HOME/pa3/PA3-Prob2.cu .

echo " "; echo "--------------------------------------------------------------"
echo "Programming Assignment 3, Problem 1"
echo "--------------------------------------------------------------"; echo " "
nvcc -O3 -arch=sm_20 PA3-Prob2.cu; ./a.out

#echo " "; echo "--------------------------------------------------------------"
#echo "Programming Assignment 3, Problem 2"
#echo "--------------------------------------------------------------"; echo " "
#nvcc -O3 -arch=sm_20 PA3-Prob2.cu; ./a.out

rm -rf $TMPDIR
