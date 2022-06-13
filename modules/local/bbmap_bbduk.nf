// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BBMAP_BBDUK {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bbmap=38.93" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bbmap%3A38.93--he522d1c_0"
    } else {
        container "quay.io/biocontainers/bbmap:bbmap:38.93--he522d1c_0"
    }

    input:
    tuple val(meta), path(reads)
    path contaminants

    output:
    tuple val(meta), path('*.qc.*fastq.gz'), emit: reads
    tuple val(meta), path('*.log')     , emit: log
    path ("versions.yml"), emit: versions


    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}.${options.suffix}" : "${meta.id}"

    def raw      = meta.single_end ? "in=${reads[0]}" : "in1=${reads[0]} in2=${reads[1]}"
    //step 1: out
    def adapter_trimmed_out  = meta.single_end ? "out=${prefix}_trimAdapter.fastq.gz" : "out1=${prefix}_trimAdapter_R1.fastq.gz out2=${prefix}_trimAdapter_R2.fastq.gz"
    
    //step2
    def q_filter_in  = meta.single_end ? "in=${prefix}_trimAdapter.fastq.gz" : "in1=${prefix}_trimAdapter_R1.fastq.gz in2=${prefix}_trimAdapter_R2.fastq.gz"
    def q_filter_out  = meta.single_end ? "out=${prefix}_qfilter.fastq.gz" : "out1=${prefix}_qfilter_R1.fastq.gz out2=${prefix}_qfilter_R2.fastq.gz"
    //step3
    def a_filter_in  = meta.single_end ? "in=${prefix}_qfilter.fastq.gz" : "in1=${prefix}_qfilter_R1.fastq.gz in2=${prefix}_qfilter_R2.fastq.gz"
    def out =  meta.single_end ? "out=${prefix}.qc.fastq.gz" : "out1=${prefix}.qc.R1.fastq.gz out2=${prefix}.qc.R2.fastq.gz"

    def maxmem = "-Xmx${task.memory.toGiga()}g"
    // https://www.protocols.io/view/illumina-fastq-filtering-gydbxs6?step=3
    """
    #Step 1: remove sequencing adapter
    #: 1)trims the last base off of 151bp reads; that base is very low quality. 
    #   Specifically, ftm=5 will trim reads so that their length is equal to zero modulo 5, 
    #and ignore reads that are already 100bp or 150bp, etc
    #2): trim off the partial adapter, artifacts, phix
    # with two step together, they do ftm=5 first and then adapter trimming

    bbduk.sh $maxmem $raw $adapter_trimmed_out $options.args threads=${task.cpus} >& ${prefix}.bbduk.step1.log

    #Step 2: quality filtering: 
    #This step removes reads and regions with low quality 
    #or match the sequencing artifacts database 

    bbduk.sh $maxmem $q_filter_in  $q_filter_out $options.args2 stats=${prefix}_stats.txt threads=${task.cpus} >& ${prefix}.bbduk.step2.log

    #Step 3: Artifact Filtering
    #This step removes reads in the short sequencing artifacts 
    #database from the fastq from STEP 2
    bbduk.sh $maxmem  $a_filter_in  $out  $options.args3 threads=${task.cpus} >& ${prefix}.bbduk.step3.log
    
     cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(bbversion.sh)
    END_VERSIONS

    """
    
}
