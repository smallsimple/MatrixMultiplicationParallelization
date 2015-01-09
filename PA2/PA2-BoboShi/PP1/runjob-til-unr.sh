#
# Batch submission script for Prog. Assignment 2
# If this script is placed in file runjob.sh,
# to submit batch job, use "qsub runjob.sh"
# Status of batch job can be checked using "qstat -u `whoami`" 
#
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=12
#PBS -N PP1
#PBS -S /bin/ksh
#PBS -j oe

TMPDIR=`mktemp -d`
cd $TMPDIR

#
# Assumes the Openmp programs are in file PA2-Prob{1,2}.c in a
# subdirectory called pa2 in your OSC home directory
# Change to match your directory and file name
#
INP=$HOME/PA2/PP1
input='PA2-Prob1-til-unr.c'
cp $INP/$input .
#cp $HOME/PA2/PA2-Prob2.c .
ogcc='R1-ikj-tilIK-unr2j-gcc.txt'
oicc='R1-ikj-tilIK-unr2j-icc.txt'

echo " "; echo "--------------------------------------------------------------"
echo "Programming Assignment 2, Problem 1, GCC"
echo "--------------------------------------------------------------"; echo " "
gcc -O3 -fopenmp $input -lm; ./a.out > $ogcc

echo " "; echo "--------------------------------------------------------------"
echo "Programming Assignment 2, Problem 1, ICC"
echo "--------------------------------------------------------------"; echo " "
icc -fast -openmp $input -lm; ./a.out >$oicc

#echo " "; echo "--------------------------------------------------------------"
#echo "Programming Assignment 2, Problem 2, GCC"
#echo "--------------------------------------------------------------"; echo " "
#gcc -O3 -fopenmp PA2-Prob2.c -lm; ./a.out

#echo " "; echo "--------------------------------------------------------------"
#echo "Programming Assignment 2, Problem 2, ICC"
#echo "--------------------------------------------------------------"; echo " "
#icc -fast -openmp PA2-Prob2.c -lm; ./a.out
cp $ogcc $oicc $INP/
rm -rf $TMPDIR
