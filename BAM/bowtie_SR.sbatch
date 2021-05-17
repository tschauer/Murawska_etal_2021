#! /bin/bash
#
#
#SBATCH -J bowtie_array   # A single job name for the array
#SBATCH -p slim18                    # Partition
#SBATCH -n 16                        # 16 cores
#SBATCH -N 1                         # one node
#SBATCH -t 0-2:00                    # Running time of 2 hours
#SBATCH --mem 20000                  # Memory request
#SBATCH -o bowtie_%A_%a.out          # Standard output
#SBATCH -e bowtie_%A_%a.err          # Standard error


FILENAME=`ls *.txt.gz  | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
FILEBASE=`echo ${FILENAME} | sed -e "s/.txt.gz//g"`

module load ngs/bowtie2/2.2.9
module load ngs/samtools/1.8


BOWTIE_INDEX="../genome/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel"

BOWTIE_OPTS="-p 16"


bowtie2 $BOWTIE_OPTS -x $BOWTIE_INDEX -U ${FILEBASE}.txt.gz > ${FILEBASE}.sam 2> ${FILEBASE}.stats


samtools view -bS -q 0 -@ 8 ${FILEBASE}.sam | samtools sort - -@ 16 | tee ${FILEBASE}.multi.bam | samtools index - ${FILEBASE}.multi.bam.bai


rm ${FILEBASE}.sam
