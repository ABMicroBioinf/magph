
/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()
//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )
include {BRACKEN} from '../modules/local/bracken'        addParams( options: modules['kraken2_bracken']) 
include {BRACKEN as BRACKEN_LONG} from '../modules/local/bracken'        addParams( options: modules['kraken2_bracken'])   
include {GET_SAMPLEIDS} from '../modules/local/get_sampleids' addParams( options: modules['csvtk_concat'])
include {METAMAPS} from '../modules/local/metamaps' addParams( options: modules['metamaps'])
include {
    CSVTK_CONCAT as CONCAT_STATS_SHORT_RAW;
    CSVTK_CONCAT as CONCAT_STATS_SHORT_QC;
    CSVTK_CONCAT as CONCAT_STATS_LONG_RAW;
    CSVTK_CONCAT as CONCAT_STATS_LONG_QC;
    CSVTK_CONCAT as CONCAT_STATS_SHORT_DEHOST;
    CSVTK_CONCAT as CONCAT_STATS_LONG_DEHOST;
} from '../modules/local/csvtk_concat' addParams( options: modules['csvtk_concat'])

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/
//
// MODULE: Installed directly from nf-core/modules
//
//include { FASTQC  } from '../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
include {KRAKEN2_KRAKEN2} from '../modules/nf-core/modules/kraken2/kraken2/main'        addParams( options: modules['kraken2'])   
include {KRAKEN2_KRAKEN2 as KRAKEN2_KRAKEN2_LONG} from '../modules/nf-core/modules/kraken2/kraken2/main'        addParams( options: modules['kraken2'])   
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main' addParams( options: [publish_files : ['_versions.yml':'']] )
include {METAPHLAN3} from '../modules/nf-core/modules/metaphlan3/main'
include {METAPHLAN3 as METAPHLAN3_LONG} from '../modules/nf-core/modules/metaphlan3/main'
//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )
include {RUN_ILLUMINA_QC} from '../subworkflows/local/qc_illumina'
include {RUN_NANOPORE_QC} from '../subworkflows/local/qc_nanopore'
include {DEHOST_ILLUMINA} from '../subworkflows/local/dehost_illumina'
include {DEHOST_NANOPORE} from '../subworkflows/local/dehost_nanopore'


/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMagph.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
//def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
def checkPathParamList = [ params.input]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

 
/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []
workflow MAGPH {

    ch_software_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    reads = INPUT_CHECK.out.reads
    short_reads = INPUT_CHECK.out.shortreads
    long_reads = INPUT_CHECK.out.longreads

    INPUT_CHECK.out.ids.collect() | GET_SAMPLEIDS

    if(!params.skip_short_reads_qc){
        RUN_ILLUMINA_QC(short_reads)
        ch_software_versions = ch_software_versions.mix(RUN_ILLUMINA_QC.out.versions)
        short_reads = RUN_ILLUMINA_QC.out.qc_reads
        CONCAT_STATS_SHORT_RAW(RUN_ILLUMINA_QC.out.raw_stats.map { cfg, stats -> stats }.collect().map { files -> tuple("short_reads_raw_seqstats", files)} )
        CONCAT_STATS_SHORT_QC(RUN_ILLUMINA_QC.out.qc_stats.map { cfg, stats -> stats }.collect().map { files -> tuple("short_reads_qc_seqstats", files)} )
        
    }
    if(!params.skip_long_reads_qc){
        RUN_NANOPORE_QC(long_reads)
        long_reads = RUN_NANOPORE_QC.out.qc_reads
        ch_software_versions = ch_software_versions.mix(RUN_NANOPORE_QC.out.versions)
        CONCAT_STATS_LONG_RAW(RUN_NANOPORE_QC.out.raw_stats.map { cfg, stats -> stats }.collect().map { files -> tuple("long_reads_raw_seqstats", files)} )
        CONCAT_STATS_LONG_QC(RUN_NANOPORE_QC.out.qc_stats.map { cfg, stats -> stats }.collect().map { files -> tuple("long_reads_qc_seqstats", files)} )
        
    }
    //remove host human reads
    //https://www.metagenomics.wiki/tools/short-read/remove-host-sequences
     if(!params.skip_remove_host){
        //BOWTIE2_BUILD (file(params.human_host))
        //DEHOST_ILLUMINA(short_reads, file(params.genomes['GRCh38']['bowtie2']))
        //use locally downloaded index
        DEHOST_ILLUMINA(short_reads, file(params.human_host_bowtie2index))
        short_reads = DEHOST_ILLUMINA.out.fastq
        ch_software_versions = ch_software_versions.mix(DEHOST_ILLUMINA.out.versions)
        CONCAT_STATS_SHORT_DEHOST(DEHOST_ILLUMINA.out.dehost_stats.map { cfg, stats -> stats }.collect().map { files -> tuple("short_reads_dehost_seqstats", files)} )
        //long reads
        //DEHOST_NANOPORE(long_reads, file(params.genomes['GRCh38']['fasta']))
        DEHOST_NANOPORE(long_reads, file(params.human_host_fasta))
        long_reads = DEHOST_NANOPORE.out.fastq
        ch_software_versions = ch_software_versions.mix(DEHOST_NANOPORE.out.versions)
        CONCAT_STATS_LONG_DEHOST(DEHOST_NANOPORE.out.dehost_stats.map { cfg, stats -> stats }.collect().map { files -> tuple("long_reads_dehost_seqstats", files)} )
    }
 
    //classify
    if(!params.skip_kraken2){
        Channel
            .value(file( "${params.kraken2_db}" ))
            .set { ch_kraken2_db_file }
        KRAKEN2_KRAKEN2 (short_reads, ch_kraken2_db_file)
        ch_software_versions = ch_software_versions.mix(KRAKEN2_KRAKEN2.out.versions)
        BRACKEN(KRAKEN2_KRAKEN2.out.txt, ch_kraken2_db_file)
        ch_software_versions = ch_software_versions.mix(BRACKEN.out.versions)

        KRAKEN2_KRAKEN2_LONG (long_reads, ch_kraken2_db_file)
        ch_software_versions = ch_software_versions.mix(KRAKEN2_KRAKEN2_LONG.out.versions)
        BRACKEN_LONG(KRAKEN2_KRAKEN2_LONG.out.txt, ch_kraken2_db_file)
        ch_software_versions = ch_software_versions.mix(BRACKEN_LONG.out.versions)
    }
    if(!params.skip_metaphlan3){
        METAPHLAN3(short_reads, file(params.metaphlan_db))
    }
    //it runs forever
    //METAMAPS(long_reads, file(params.metamaps_db))

    /* CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_software_versions.unique().collectFile()
        //ch_software_versions.collectFile()
    ) */
} 

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
