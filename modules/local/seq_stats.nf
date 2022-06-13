
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'
params.options = [:]
options        = initOptions(params.options)

process SEQ_STATS{

    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
    
    input:
    tuple val(meta), path(fqchk_out)

    output:
    tuple val(meta), path('*_seqstats.txt'), emit: stats

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.suffix ? "${meta.id}.${options.suffix}" : "${meta.id}"
        def seqtype = meta.single_end ? "single" : "paired"
        if( options.args == "long"){
            seqtype = "single"
        }

        """
        parse_seqtk_fqchk.pl  -t  ${seqtype} -i ${fqchk_out} -s ${meta.id} > ${prefix}_seqstats.txt
        """
   
}
