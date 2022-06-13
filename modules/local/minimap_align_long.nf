// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'


params.options = [:]
options        = initOptions(params.options)

process MINIMAP2_ALIGN_LONG {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
    

    conda (params.enable_conda ? 'bioconda::minimap2=2.22' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/minimap2:2.21--h5bf99c6_0"
    } else {
        container "quay.io/biocontainers/minimap2:2.21--h5bf99c6_0"
    }

    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(reference)

    output:
    tuple val(meta), path("*.paf"), emit: paf
    path "versions.yml" , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    
    """
    minimap2 \\
        $options.args \\
        -t $task.cpus \\
        $reference \\
        $reads \\
        > ${prefix}.paf

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
