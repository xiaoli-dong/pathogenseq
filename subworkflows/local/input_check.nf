//
// Check input samplesheet and get read channels
//

include { SAMPLESHEETCHECK } from '../../modules/local/samplesheetcheck'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEETCHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_read_channels(it) }
        .set { reads }

    //reads.view()

    reads.map { 
        meta, reads, long_fastq, fast5 -> [ meta, reads ] }
        .filter{ meta, reads -> reads[0] != 'NA' && reads[1] != 'NA' }
        .set { shortreads }
    
    //shortreads.view() 
    

    reads.map {
        meta, reads, long_fastq, fast5 -> [ meta, long_fastq ] }
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


    //emit:
    //reads                                     // channel: [ val(meta), [ reads ] ]
    //versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
    
    emit:
    reads      // channel: [ val(meta), [ reads ], long_fastq, fast5 ]
    shortreads // channel: [ val(meta), [ reads ] ]
    longreads  // channel: [ val(meta), long_fastq ]
    fast5      // channel: [ val(meta), fast5 ]
    ids
    versions = SAMPLESHEETCHECK.out.versions // channel: [ versions.yml ]
}
// Function to get list of [ meta, [ fastq_1, fastq_2 ], long_fastq, fast5 ]
def create_read_channels(LinkedHashMap row) {
    
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = !(row.fastq_1 == 'NA') && !(row.fastq_2 == 'NA') ? false : true
    meta.basecaller_mode = row.basecaller_mode == null ? 'NA' : row.basecaller_mode
    meta.genome_size  = row.genomesize == null ? 'NA' : row.genomesize
    
    if( !(row.fastq_1 == 'NA') )
    def array = []
    // check short reads
    if ( !(row.fastq_1 == 'NA') ) {
        if ( !file(row.fastq_1).exists() ) {
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
        }
        
        if(file(row.fastq_1).size() == 0){
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file is empty!\n${row.fastq_1}" 
        }
        
        fastq_1 = file(row.fastq_1)
    } else { fastq_1 = 'NA' }
    if ( !(row.fastq_2 == 'NA') ) {
        if ( !file(row.fastq_2).exists() ) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        
        if(file(row.fastq_2).size() == 0){
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file is empty!\n${row.fastq_2}" 
        }
        fastq_2 = file(row.fastq_2)
    } else { fastq_2 = 'NA' }

    // check long_fastq
    if ( !(row.long_fastq == 'NA') ) {
        if ( !file(row.long_fastq).exists() ) {
            exit 1, "ERROR: Please check input samplesheet -> Long FastQ file does not exist!\n${row.long_fastq}"
        }
        
        if(file(row.long_fastq).size() == 0){
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file is empty!\n${row.long_fastq}" 
        }
        long_fastq = file(row.long_fastq)
    } else { long_fastq = 'NA' }

    // check fast5
    if ( !(row.fast5 == 'NA') ) {
        if ( !file(row.fast5).exists() ) {
            exit 1, "ERROR: Please check input samplesheet -> fast5 file does not exist!\n${row.fast5}"
        }
       
        if(file(row.fast5).size() == 0){
            exit 1, "ERROR: Please check input samplesheet -> fast5 file is empty!\n${row.fast5}" 
        }

        fast5 = file(row.fast5)
    } else { fast5 = 'NA' }

    // check basecaller_mode
    def basecaller_modes = ["fast", "hac", "sup"]
    if(! meta.basecaller_mode.toLowerCase() in basecaller_modes){
        meta.basecaller_mode = 'NA'
    }
    
    // prepare output // currently does not allow single end data!
    if ( meta.single_end ) {
        array = [ meta, fastq_1 , long_fastq, fast5 ]
    } else {
        array = [ meta, [ fastq_1, fastq_2 ], long_fastq, fast5 ]
    } 
    return array 
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel_org(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}
