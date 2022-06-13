#!/usr/bin/env nextflow

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//nanopore
include {NANOPLOT as NANOPLOT_RAW} from '../../modules/nf-core/modules/nanoplot/main' addParams( options: modules['nanoplot_raw']) 
include {NANOPLOT as NANOPLOT_QC}  from '../../modules/nf-core/modules/nanoplot/main' addParams( options: modules['nanoplot_qc']) 
include {PORECHOP} from '../../modules/local//porechop'                              addParams( options: modules['porechop']) 
include {FILTLONG} from '../../modules/local//filtlong'                              addParams( options: modules['filtlong']) 

include {SEQTK_FQCHK as SEQTK_FQCHK_RAW} from '../../modules/local/seqtk_fqchk' addParams( options: modules['seqtk_fqchk_longreads_raw']) 
include {SEQ_STATS as SEQ_STATS_RAW} from '../../modules/local/seq_stats' addParams( options: modules['seq_stats_longreads_raw']) 
include {SEQTK_FQCHK as SEQTK_FQCHK_QC} from '../../modules/local/seqtk_fqchk' addParams( options: modules['porechop_qc_seqtk_fqchk']) 
include {SEQ_STATS as SEQ_STATS_QC} from '../../modules/local/seq_stats' addParams( options: modules['porechop_qc_seq_stats']) 


workflow RUN_NANOPORE_QC {

    take:
        reads
    main:

        ch_versions = Channel.empty()

        /* NANOPLOT_RAW(reads)
        ch_versions = ch_versions.mix(NANOPLOT_RAW.out.versions.first()) */
        
        SEQTK_FQCHK_RAW(reads)
        ch_versions = ch_versions.mix(SEQTK_FQCHK_RAW.out.versions.first())
        SEQ_STATS_RAW(SEQTK_FQCHK_RAW.out.stats)

        // NANOPORESTATS_RAW(reads)
        // raw_stats = NANOPORESTATS_RAW.out.stats
        //ch_versions = ch_versions.mix(NANOPORESTATS_RAW.out.versions.first())
        PORECHOP(reads)

        porechop_reads = PORECHOP.out.reads
        FILTLONG(porechop_reads)
        
        
        qc_reads = FILTLONG.out.reads
        ch_versions = ch_versions.mix(PORECHOP.out.versions.first())
        SEQTK_FQCHK_QC(qc_reads)
        SEQ_STATS_QC(SEQTK_FQCHK_QC.out.stats)

        // NANOPORESTATS_QC(reads)
        // qc_stats = NANOPORESTATS_QC.out.stats


        //NANOPLOT_QC(qc_reads)


    emit:
        raw_stats = SEQ_STATS_RAW.out.stats
        qc_reads
        qc_stats = SEQ_STATS_QC.out.stats
        versions = ch_versions

}
