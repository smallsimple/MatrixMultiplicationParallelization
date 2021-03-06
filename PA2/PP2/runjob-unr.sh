#
# Batch submission script for Prog. Assignment 2
# If this script is placed in file runjob.sh,
# to submit batch job, use "qsub runjob.sh"
# Status of batch job can be checked using "qstat -u `whoami`" 
#
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=12
#PBS -N PP2
#PBS -S /bin/ksh
#PBS -j oe

TMPDIR=`mktemp -d`
cd $TMPDIR

#
# Assumes the Openmp programs are in file PA2-Prob{1,2}.c in a
# subdirectory called pa2 in your OSC home directory
# Change to match your directory and file name
#
INP=$HOME/PA2/PP2
input='PA2-Prob2-unr.c'
#cp $HOME/PA2/PA2-Prob1.c .
cp $INP/$input .
output='R2-klij-unr-4j.txt'

#echo " "; echo "--------------------------------------------------------------"
#echo "Programming Assignment 2, Problem 1, GCC"
#echo "--------------------------------------------------------------"; echo " "
#gcc -O3 -fopenmp PA2-Prob1.c -lm; ./a.out

#echo " "; echo "--------------------------------------------------------------"
#echo "Programming Assignment 2, Problem 1, ICC"
#echo "--------------------------------------------------------------"; echo " "
#icc -fast -openmp PA2-Prob1.c -lm; ./a.out

echo " "; echo "--------------------------------------------------------------"
echo "Programming Assignment 2, Problem 2, GCC">$output
echo "--------------------------------------------------------------"; echo " "
gcc -O3 -fopenmp $input -lm; ./a.out >> $output

echo " "; echo "--------------------------------------------------------------"
echo " ">>$output;echo " ">>$output;
echo "Programming Assignment 2, Problem 2, ICC">>$output
echo "--------------------------------------------------------------"; echo " "
icc -fast -openmp $input -lm; ./a.out>> $output

cp $output $INP
rm -rf $TMPDIR
