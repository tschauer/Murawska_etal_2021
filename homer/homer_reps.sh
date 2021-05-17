#!/bin/sh
# homer.sh


### makeTagDirectory ###

makeTagDirectory H3K9me2_wt_1.dir -mapq 1 ../BAM/H3K9me2_wt1_AGTCAA_1.multi.bam
makeTagDirectory H3K9me2_wt_2.dir -mapq 1 ../BAM/H3K9me2_wt2_AGTTCC_1.multi.bam

makeTagDirectory  IN_wt_1.dir -mapq 1 ../BAM/IN_wt1_ATCACG_1.multi.bam
makeTagDirectory  IN_wt_2.dir -mapq 1 ../BAM/IN_wt2_CGATGT_1.multi.bam


# peak finding against input

findPeaks  H3K9me2_wt_1.dir/ -i  IN_wt_1.dir/ -style histone -F 2 -o H3K9me2_wt_1.histone.F2.txt
findPeaks  H3K9me2_wt_2.dir/ -i  IN_wt_2.dir/ -style histone -F 2 -o H3K9me2_wt_2.histone.F2.txt

pos2bed.pl H3K9me2_wt_1.histone.F2.txt > H3K9me2_wt_1.histone.F2.bed
pos2bed.pl H3K9me2_wt_2.histone.F2.txt > H3K9me2_wt_2.histone.F2.bed
