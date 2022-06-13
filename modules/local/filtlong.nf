
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options    = initOptions(params.options)
process FILTLONG {
    tag "$meta.id"
    label 'process_medium'
     publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::filtlong=0.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'https://depot.galaxyproject.org/singularity/filtlong:0.2.1--h9a82719_0' :'quay.io/biocontainers/filtlong:0.2.1--h9a82719_0' }"

    input:
    tuple val(meta),path(longreads)

    output:
    tuple val(meta), path("*_filtlong.fastq.gz"), emit: reads
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
   
    """
    filtlong \\
        $options.args \\
        $longreads \\
        | gzip -n > ${prefix}_filtlong.fastq.gz
        
    cat <<-END_VERSIONS > versions.yml

    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( filtlong --version | sed -e "s/Filtlong v//g" )

    END_VERSIONS
    
    """
}