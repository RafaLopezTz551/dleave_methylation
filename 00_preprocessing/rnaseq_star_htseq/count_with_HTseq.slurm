#!/bin/bash

#SBATCH --job-name=singlecpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

module load anaconda3/2025.06 samtools/1.22.1
source activate htseq-2.0.9
# Your script goes here
#fungicida=($(find . -maxdepth 1 -name "l*.Aligned.out.bam"))
tentacles=($(find . -maxdepth 1 -name "*_namesorted.bam"))
for bam_file in ${tentacles[@]}
do
sorted_name=$(basename -s "_namesorted.bam" $bam_file)
#samtools sort -n -@ 8 $bam_file -O BAM -o ${sorted_name}_namesorted.bam
result_file="eye_${sorted_name}_htseq_gene_counts.txt"
#For paired end they must be sorted by name
#htseq-count -f bam -r name --stranded=reverse ${sorted_name}_namesorted.bam evian.gff3 > counts_HTseq_EviAnn/$result_file
#htseq-count -f bam --stranded=reverse $bam_file evian.gff3 > counts_HTseq_EviAnn/$result_file
htseq-count -f bam -r name --stranded=reverse $bam_file reference_coding_and_lncRNAs.gtf > counts_HTseq/$result_file
done

