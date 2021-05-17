#! /bin/bash



### run bowtie array ###


FILES=($(ls -1 *.txt.gz))

# get size of array
NUMFASTQ=${#FILES[@]}

# now submit to SLURM
if [ $NUMFASTQ -ge 0 ]; then
	sbatch --array=1-$NUMFASTQ bowtie_SR.sbatch
fi
