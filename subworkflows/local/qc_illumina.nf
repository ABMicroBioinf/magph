#!/usr/bin/env nextflow

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include {FASTQC as FASTQC_RAW} from     '../../modules/nf-core/modules/fastqc/main'         addParams( options: modules['fastqc_raw'] )
include {FASTQC as FASTQC_QC} from      '../../modules/nf-core/modules/fastqc/main'         addParams( options: modules['fastqc_qc']) 
include { BBMAP_BBDUK as BBDUK } from '../../modules/local/bbmap_bbduk'                       addParams( options: modules['bbmap_bbduk']) 

include {SEQTK_FQCHK as SEQTK_FQCHK_RAW} from '../../modules/local/seqtk_fqchk' addParams( options: modules['seqtk_fqchk_shortreads_raw']) 
include {SEQTK_FQCHK as SEQTK_FQCHK_QC} from '../../modules/local/seqtk_fqchk' addParams( options: modules['bbmap_bbduk_qc_seqtk_fqchk']) 
include {SEQ_STATS as SEQ_STATS_RAW} from '../../modules/local/seq_stats' addParams( options: modules['seq_stats_shortreads_raw']) 
include {SEQ_STATS as SEQ_STATS_QC} from '../../modules/local/seq_stats' addParams( options: modules['bbmap_bbduk_qc_stats']) 



workflow RUN_ILLUMINA_QC {   

    take:
        reads
    main:
        ch_versions = Channel.empty()
        FASTQC_RAW(reads)
        ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())
        SEQTK_FQCHK_RAW(reads)
        ch_versions = ch_versions.mix(SEQTK_FQCHK_RAW.out.versions.first())
        SEQ_STATS_RAW(SEQTK_FQCHK_RAW.out.stats)
        
        
        BBDUK(reads, [])
        ch_versions = ch_versions.mix(BBDUK.out.versions.first())
        qc_reads = BBDUK.out.reads
        FASTQC_QC(qc_reads)
        SEQTK_FQCHK_QC(qc_reads)
        ch_versions = ch_versions.mix(SEQTK_FQCHK_QC.out.versions.first())
        SEQ_STATS_QC(SEQTK_FQCHK_QC.out.stats)
        
    emit:
        qc_reads
        versions = ch_versions
        raw_stats = SEQ_STATS_RAW.out.stats
        qc_stats = SEQ_STATS_QC.out.stats
        
}

