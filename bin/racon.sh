 #!/bin/env bash

genome=$1
reads=$2
threads=$3
extra_params = $4

#Set input and parameters
round=4
input=${genome}
for ((i=1; i<=${round};i++)); do
#step 1:
    minimap2 -x map-ont -t ${threads}  ${input} ${reads} > racon_${i}.paf
    racon $extra_params -t ${threads}  $reads racon_${i}.paf ${input} > racon_round_${i}.fasta
    input=racon_round_${i}.fasta;
done;

cp racon_round_${round}.fasta genome.racon.fasta
