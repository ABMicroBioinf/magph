// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process PORECHOP {
    tag "$meta.id"

    label 'process_medium'
    
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::porechop=0.2.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/porechop%3A0.2.4--py36h2ad2d48_3"
    } else {
        container "quay.io/biocontainers/porechop:0.2.4--py36h2ad2d48_3"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_porechop.fastq")  , emit: reads
    path 'versions.yml', emit: versions

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    porechop $options.args -i ${reads} -t ${task.cpus} -o ${prefix}_porechop.fastq

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(porechop --version 2>&1)
    END_VERSIONS


    """
}