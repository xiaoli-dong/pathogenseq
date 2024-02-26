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
        meta, reads, long_fastq -> [ meta, reads ] }//.view()
        .filter{ meta, reads -> !reads[0].equals('NA') && !reads[1].equals('NA') }
        .set { shortreads }
    
    //shortreads.view() 
    

    reads.map {
        meta, reads, long_fastq -> [ meta, long_fastq ] }
        .filter{ meta, long_fastq -> !long_fastq.equals('NA') }
        .set { longreads }
       
    //longreads.view()
    
    
    reads
        .map { meta, reads, long_fastq -> meta.id }
        .set {ids}

    //ids.view()


    //emit:
    //reads                                     // channel: [ val(meta), [ reads ] ]
    //versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
    
    emit:
    reads      // channel: [ val(meta), [ reads ], long_fastq ]
    shortreads // channel: [ val(meta), [ reads ] ]
    longreads  // channel: [ val(meta), long_fastq ]
    ids
    versions = SAMPLESHEETCHECK.out.versions // channel: [ versions.yml ]
}
// Function to get list of [ meta, [ fastq_1, fastq_2 ], long_fastq ]
def create_read_channels(LinkedHashMap row) {
    
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = !(row.fastq_1 == 'NA') && !(row.fastq_2 == 'NA') ? false : true
    meta.basecaller_mode = row.basecaller_mode == null ? 'NA' : row.basecaller_mode
    
    
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

   /*  // check basecaller_mode
    //def basecaller_modes = ["fast", "hac", "sup"]
    def medaka_current_models = [
       
        //r1041 e82 (kit14) consensus
        'r1041_e82_400bps_hac_v4.2.0',
        'r1041_e82_400bps_sup_v4.2.0',
    ]
    medaka_archived_models = [
        //r9 consensus
        'r941_sup_plant_g610',
        'r941_min_fast_g507', 'r941_prom_fast_g507',
        'r941_min_fast_g303', 'r941_min_high_g303', 'r941_min_high_g330',
        'r941_prom_fast_g303', 'r941_prom_high_g303', 'r941_prom_high_g330',
        'r941_min_high_g344', 'r941_min_high_g351', 'r941_min_high_g360',
        'r941_prom_high_g344', 'r941_prom_high_g360', 'r941_prom_high_g4011',
        //r10 consensus
        'r10_min_high_g303', 'r10_min_high_g340',
        'r103_min_high_g345', 'r103_min_high_g360', 'r103_prom_high_g360',
        'r103_fast_g507', 'r103_hac_g507', 'r103_sup_g507',
        //r104 e81 consensus
        'r104_e81_fast_g5015', 'r104_e81_sup_g5015', 'r104_e81_hac_g5015',
        'r104_e81_sup_g610',
       
        //r1041 e82 consensus
        'r1041_e82_400bps_hac_g615',  'r1041_e82_400bps_fast_g615',
        'r1041_e82_400bps_fast_g632', 'r1041_e82_260bps_fast_g632',
        'r1041_e82_400bps_hac_g632', 'r1041_e82_400bps_sup_g615',
        'r1041_e82_260bps_hac_g632', 'r1041_e82_260bps_sup_g632',
        'r1041_e82_400bps_hac_v4.0.0', 'r1041_e82_400bps_sup_v4.0.0',
        'r1041_e82_260bps_hac_v4.0.0', 'r1041_e82_260bps_sup_v4.0.0',
        'r1041_e82_260bps_hac_v4.1.0', 'r1041_e82_260bps_sup_v4.1.0',
        'r1041_e82_400bps_hac_v4.1.0', 'r1041_e82_400bps_sup_v4.1.0',
        
        //rle consensus
        'r941_min_high_g340_rle',
        //r9 consensus
        'r941_min_hac_g507', 'r941_min_sup_g507',
        'r941_prom_hac_g507', 'r941_prom_sup_g507',
        
        //r941 e81 consensus
        'r941_e81_fast_g514', 'r941_e81_hac_g514', 'r941_e81_sup_g514'
       
    ]
    def medaka_allowed_models = medaka_current_models + medaka_archived_models
    //when basecaller mode missing, use medaka default mode
    if(! meta.basecaller_mode.toLowerCase() in medaka_allowed_models){
        meta.basecaller_mode = 'r1041_e82_400bps_sup_v4.2.0'
    }
     */
    // prepare output // currently does not allow single end data!
    if ( meta.single_end ) {
        array = [ meta, [fastq_1, fastq_2] , long_fastq]
    } else {
        array = [ meta, [ fastq_1, fastq_2 ], long_fastq]
    } 
    return array 
}
