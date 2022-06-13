
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SEQTK_FQCHK{

    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
    
    conda (params.enable_conda ? "bioconda::seqtk=1.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3"
    } else {
        container "quay.io/biocontainers/seqtk:1.3--h5bf99c6_3"
    }


    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*_seqtk.txt'), emit: stats
    path ("versions.yml"), emit: versions


    script:
    
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}.${options.suffix}" : "${meta.id}"
    def input      = meta.single_end ? "${reads[0]}" : "${reads[0]} ${reads[1]}"
    
    if( options.args == "long"){
            input = reads[0]
        }
    
    if(reads[0] =~ /\.gz$/)
      """
        zcat ${input} | seqtk fqchk -q0  - > ${prefix}_seqtk.txt

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(seqtk 2>&1 | grep Version | sed 's/^Version: //')
        END_VERSIONS
        """

    else 

    """
        cat ${input} | seqtk fqchk -q0  - > ${prefix}_seqtk.txt
        
        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(seqtk 2>&1 | grep Version | sed 's/^Version: //')
        END_VERSIONS
    """
    
}
