// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process METAMAPS {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::metamaps=0.1.98102e9' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/metamaps%3A0.1.98102e9--h176a8bc_0'
    } else {
        container 'quay.io/biocontainers/metamaps:0.1.98102e9--h176a8bc_0'
    }

    input:
    tuple val(meta), path(reads)
    path  db

    output:
    tuple val(meta), path('*WIMP')  , emit: WIMP
    tuple val(meta), path('*reads2Taxon'), emit: reads2Taxon
    tuple val(meta), path('*krona')   , emit: krona
    tuple val(meta), path('*contigCoverage')   , emit: contigCoverage
    tuple val(meta), path('*lengthAndIdentitiesPerMappingUnit')   , emit: lengthAndIdentitiesPerMappingUnit
    tuple val(meta), path('*EM')   , emit: EM
    tuple val(meta), path('*evidenceUnknownSpecies')   , emit: evidenceUnknownSpecies
    path "versions.yml"                    , emit: versions

    script:
    def prefix       = options.suffix  ? "${meta.id}${options.suffix}"  : "${meta.id}"
    
    """
    metamaps mapDirectly --all -t $task.cpus -r ${db}/DB.fa -q $reads -o classification_results --maxmemory 12
    metamaps classify -t $task.cpus --mappings classification_results --DB ${db}

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(metamaps 2>&1 |head -2 | sed '1d' | sed 's/^MetaMaps v //;' )
        
    END_VERSIONS
    """
}
