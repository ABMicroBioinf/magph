#!/usr/bin/bash

genome=$1
read1=$2
read2=$3
threads=$4
maxmem=$5
$extra_params = $6

pilon=/data/software/miniconda3/envs/pilon
#samtools, minimap2 are required

#Set input and parameters
round=4
input=${genome}
for ((i=1; i<=${round};i++)); do
#step 1:
        
    minimap2 -ax sr -t ${threads}  ${input} ${read1} ${read2} | samtools view --threads ${threads} -F 0x4 -b - | samtools sort - -m 2g --threads ${threads} -o sgs.sort.bam;
    #index bam and genome files
    samtools index -@ $task.cpus sgs.sort.bam;
    samtools faidx ${input};
    #polish genome file
    java ${maxmem} -jar  ${pilon}/share/pilon-1.24-0/pilon.jar --genome ${input} --frags sgs.sort.bam --output pilon_round_${i} $extra_params
    input=pilon_round_${i}.fasta;
done;

cp pilon_round_${round}.fasta genome.pilon.fasta
#Finally polished genome file: genome.pilon.fa