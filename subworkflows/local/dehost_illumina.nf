#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
def modules = params.modules.clone()
include { BOWTIE2_ALIGN     } from '../../modules/nf-core/modules//bowtie2/align/main' addParams( options: modules['bowtie2_align']) 
include { SAMTOOLS_SORT      } from '../../modules/nf-core/modules/samtools/sort/main'  addParams( options: modules['samtools_sort']) 
include { SAMTOOLS_INDEX     } from '../../modules/nf-core/modules/samtools/index/main' addParams( options: modules['samtools_index'] )
include { SAMTOOLS_VIEW      } from '../../modules/nf-core/modules/samtools/view/main'  addParams( options: modules['samtools_view_dehost_short']) 
include { SAMTOOLS_FASTQ    } from '../../modules/nf-core/modules//samtools/fastq/main' 

include {SEQTK_FQCHK as SEQTK_FQCHK_DEHOST} from '../../modules/local/seqtk_fqchk' addParams( options: modules['seqtk_fqchk_shortreads_dehost'])
include {SEQ_STATS as SEQ_STATS_DEHOST} from '../../modules/local/seq_stats' addParams( options: params.modules['seq_stats_shortreads_dehost'])

workflow DEHOST_ILLUMINA {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    index // channel: /path/to/bowtie2/index/

    main:
    ch_versions = Channel.empty()
    print "******************************************xiaoli"
    
    reads.view()
    
    //
    // Map reads with Bowtie2 to host genome and return dehosted bam
    //
    // do not save unmapped reads, only covert output sam file to bam
    BOWTIE2_ALIGN ( reads, index, false, false )
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

    //filter out the mapped human reads
    SAMTOOLS_VIEW(BOWTIE2_ALIGN.out.bam)    //

    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    SAMTOOLS_SORT (  SAMTOOLS_VIEW.out.bam )
    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    SAMTOOLS_FASTQ(SAMTOOLS_SORT.out.bam)
    
    SEQTK_FQCHK_DEHOST(SAMTOOLS_FASTQ.out.fastq)
    ch_versions = ch_versions.mix(SEQTK_FQCHK_DEHOST.out.versions.first())
    SEQ_STATS_DEHOST(SEQTK_FQCHK_DEHOST.out.stats)

    emit:
    dehost_stats = SEQ_STATS_DEHOST.out.stats
    fastq             = SAMTOOLS_FASTQ.out.fastq
    versions       = ch_versions                      // channel: [ versions.yml ]
}