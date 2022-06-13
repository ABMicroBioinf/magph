#!/usr/bin/bash

genome=$1
read1=$2
read2=$3
threads=$4
extra_params=$5

NextPolish=/data/software/NextPolish
#samtools, minimap2 are required

#Set input and parameters
round=2
input=${genome}
for ((i=1; i<=${round};i++)); do

#two different algorithm
#step 1:
        minimap2 -ax sr -t ${threads}  ${input} ${read1} ${read2} | samtools view --threads ${threads} -F 0x4 -b - | samtools sort - -m 2g --threads ${threads} -o sgs.sort.bam;
        #index bam and genome files
        samtools index -@ 20 sgs.sort.bam;
        samtools faidx ${input};
        #polish genome file
        python $NextPolish/lib/nextpolish1.py -g ${input} -t 1 -p ${threads} -s sgs.sort.bam ${extra_params}> genome.polishtemp.fasta;
        input=genome.polishtemp.fasta;
#step2:
        minimap2 -ax sr -t ${threads}  ${input} ${read1} ${read2} | samtools view --threads ${threads} -F 0x4 -b - | samtools sort - -m 2g --threads ${threads} -o sgs.sort.bam;
        #index bam and genome files
        samtools index -@ 20 sgs.sort.bam;
        samtools faidx ${input};
        #polish genome file
        python $NextPolish/lib/nextpolish1.py -g ${input} -t 2 -p ${threads} -s sgs.sort.bam ${extra_params}> genome.nextpolish.fasta;
        input=genome.nextpolish.fasta;
done;
#Finally polished genome file: genome.nextpolish.fasta