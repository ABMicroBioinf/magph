// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GET_SAMPLEIDS{
     publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

   
    input:
    val ids

    output:
    path("sids.txt"), emit: sample_ids

    script:
    def input = ids.join("\n")
    """
    echo "${input}" > sids.txt

    """
}

