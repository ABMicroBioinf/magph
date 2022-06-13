// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CSVTK_CONCAT {

    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::csvtk=v0.23.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/csvtk%3A0.23.0--h9ee0642_0"
    } else {
        container "quay.io/biocontainers/csvtk:0.23.0--h9ee0642_0"
    }

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
    
    
    input:
    tuple val(output_name), val(input)

    output:
    path("*.tsv"), emit: collated
    path  'versions.yml',             emit: versions

    script:
    def input_files = input.join(' ')
    """
    csvtk concat -t -T $input_files > ${output_name}.tsv

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(csvtk version)
    END_VERSIONS

    """
        
}
