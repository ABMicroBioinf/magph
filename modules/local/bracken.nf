// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BRACKEN {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bracken=2.6.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bracken%3A2.6.1--py39h7cff6ad_2"
    } else {
        container "quay.io/biocontainers/bracken:2.6.1--py39h7cff6ad_2"
    }

    input:
    tuple val(meta), path(kraken2_report)
    path kraken2_db

    output:
    tuple val(meta), path('*_bracken.output.txt'), emit: output
    tuple val(meta), path('*_bracken.outreport.txt'), emit: outreport
    path ("versions.yml"), emit: versions

    script:
    def prefix  = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
   
    """
    bracken \\
        $options.args \\
        -t $task.cpus \\
        -d $kraken2_db \\
        -i ${kraken2_report} \\
        -o ${prefix}_bracken.output.txt \\
        -w ${prefix}_bracken.outreport.txt 
   
    printf "BRACKEN:\n  bracken: 2.6.1\n" > versions.yml
    """
}
