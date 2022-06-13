//
// Check input samplesheet and get read channels
//

params.options = [:]

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:

    println samplesheet
    
    Channel
        .fromPath( samplesheet )
        .ifEmpty {exit 1, log.info "Cannot find path file ${tsvFile}"}
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channels(it) }
        .set { reads } // channel: [ val(meta), [ reads ] ]

    println(reads)
        
    // reconfigure channels
    reads
        .map { meta, reads, long_fastq, fast5 -> [ meta, reads ] }
        .filter{ meta, reads -> reads != 'NA' }
        .filter{ meta, reads -> reads[0] != 'NA' && reads[1] != 'NA' }
        .set { shortreads }
    
    //shortreads.view() 
    

    reads
        .map { meta, reads, long_fastq, fast5 -> [ meta, long_fastq ] }
        .filter{ meta, long_fastq -> long_fastq != 'NA' }
        .set { longreads }
       
    //longreads.view() 
    
    reads
        .map { meta, reads, long_fastq, fast5 -> [ meta, fast5 ] }
        .filter{ meta, fast5 -> fast5 != 'NA' }
        .set { fast5 }
       

    //fast5.view() 

    
    reads
        .map { meta, reads, long_fastq, fast5 -> meta.id }
        
        .set {ids}

    //ids.view()
    

    emit:
    reads      // channel: [ val(meta), [ reads ], long_fastq, fast5 ]
    shortreads // channel: [ val(meta), [ reads ] ]
    longreads  // channel: [ val(meta), long_fastq ]
    fast5      // channel: [ val(meta), fast5 ]
    ids
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ], long_fastq, fast5 ]
def create_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    
    meta.id           = row.sample
    meta.single_end   = !(row.fastq_1 == 'NA') && !(row.fastq_2 == 'NA') ? false : true

    meta.genome_size  = row.GenomeSize == null ? 'NA' : row.GenomeSize
    
    if( !(row.fastq_1 == 'NA') )
    def array = []
    // check short reads
    if ( !(row.fastq_1 == 'NA') ) {
        if ( !file(row.fastq_1).exists() ) {
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
        }
        fastq_1 = file(row.fastq_1)
    } else { fastq_1 = 'NA' }
    if ( !(row.fastq_2 == 'NA') ) {
        if ( !file(row.fastq_2).exists() ) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_2 = file(row.fastq_2)
    } else { fastq_2 = 'NA' }

    // check long_fastq
    if ( !(row.long_fastq == 'NA') ) {
        if ( !file(row.long_fastq).exists() ) {
            exit 1, "ERROR: Please check input samplesheet -> Long FastQ file does not exist!\n${row.long_fastq}"
        }
        long_fastq = file(row.long_fastq)
    } else { long_fastq = 'NA' }

    // check long_fastq
    if ( !(row.Fast5 == 'NA') ) {
        if ( !file(row.Fast5).exists() ) {
            exit 1, "ERROR: Please check input samplesheet -> Fast5 file does not exist!\n${row.Fast5}"
        }
        fast5 = file(row.Fast5)
    } else { fast5 = 'NA' }

    // prepare output // currently does not allow single end data!
    if ( meta.single_end ) {
        array = [ meta, fastq_1 , long_fastq, fast5 ]
    } else {
        array = [ meta, [ fastq_1, fastq_2 ], long_fastq, fast5 ]
    } 
    return array
}

