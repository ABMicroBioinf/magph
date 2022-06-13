#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
def modules = params.modules.clone()

include {MINIMAP2_ALIGN} from '../../modules/local/minimap2_align' 
include { SAMTOOLS_VIEW      } from '../../modules/nf-core/modules/samtools/view/main'  addParams( options: modules['samtools_view_dehost_long']) 
include { SAMTOOLS_SORT      } from '../../modules/nf-core/modules/samtools/sort/main'  addParams( options: modules['samtools_sort_long']) 
include { SAMTOOLS_INDEX     } from '../../modules/nf-core/modules/samtools/index/main' addParams( options: modules['samtools_index_long'] )
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_LONG    } from '../../modules/nf-core/modules//samtools/fastq/main' 
include {SEQTK_FQCHK as SEQTK_FQCHK_DEHOST} from '../../modules/local/seqtk_fqchk' addParams( options: modules['seqtk_fqchk_longreads_dehost'])
include {SEQ_STATS as SEQ_STATS_DEHOST} from '../../modules/local/seq_stats' addParams( options: params.modules['seq_stats_longreads_dehost'])

workflow DEHOST_NANOPORE {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    reference 
    main:
    ch_versions = Channel.empty()

    //
    // Map reads with Bowtie2 to host genome and return dehosted bam
    //
    MINIMAP2_ALIGN ( reads, reference)
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    //filter out the mapped human reads
    SAMTOOLS_VIEW(MINIMAP2_ALIGN.out.sam)    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    SAMTOOLS_SORT ( SAMTOOLS_VIEW.out.bam )
    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    SAMTOOLS_FASTQ_LONG(SAMTOOLS_SORT.out.bam)
    
    SEQTK_FQCHK_DEHOST(SAMTOOLS_FASTQ_LONG.out.fastq)
    ch_versions = ch_versions.mix(SEQTK_FQCHK_DEHOST.out.versions.first())
    SEQ_STATS_DEHOST(SEQTK_FQCHK_DEHOST.out.stats)

    emit:
    dehost_stats = SEQ_STATS_DEHOST.out.stats
    fastq             = SAMTOOLS_FASTQ_LONG.out.fastq
    versions       = ch_versions                      // channel: [ versions.yml ]
}