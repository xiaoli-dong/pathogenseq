//
// prepare databases
//

include { AMRFINDERPLUS_UPDATE          } from '../../modules/nf-core/amrfinderplus/update/main'

workflow PREPARE_REFERENCES {
    main:

    ch_versions = Channel.empty()

    
    //
    // Prepare kraken db
    //
    ch_kraken2_db = Channel.empty()
    if (!params.skip_illumina_kraken2) {
        if (params.kraken2_db) {
            ch_kraken2_db = Channel.value(file(params.kraken2_db))
        }
        else {
            log.error "Please specify a valid  Kraken2 database."
            System.exit(1)
        }
    }

    //
    // Prepare kraken db
    //
    ch_bakta_db = Channel.empty()
    if (!params.skip_bakta) {
        if (params.bakta_db) {
            ch_bakta_db = Channel.value(file(params.bakta_db))
        }
        else {
            log.error "Please specify a valid  bakta database."
            System.exit(1)
        }
    }

    //
    // Prepare kraken db
    //
    ch_checkm2_db = Channel.empty()
    if (!params.skip_checkm2) {
        if (params.checkm2_db) {
            ch_checkm2_db = Channel.value(file(params.checkm2_db))
        }
        else {
            log.error "Please specify a valid checkm2 database."
            System.exit(1)
        }
    }

    //
    // Prepare amrfinderplus database
    //
    ch_amrfinderplus_db = Channel.empty()
    if (!params.skip_amr) {
        if (params.amrfinderplus_db) {
            ch_checkm2_db = Channel.value(file(params.checkm2_db))
        }
        else {
            AMRFINDERPLUS_UPDATE()
            ch_amrfinderplus_db = AMRFINDERPLUS_UPDATE.out.db
            ch_versions   = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions)
        }
    }

    
    emit:
        ch_kraken2_db                    // path: kraken2_db/
        ch_checkm2_db   //path to checkm2 db
        ch_bakta_db     // path to bakta database 
        ch_amrfinderplus_db // path: path to amrfinderplus_db
        versions             = ch_versions             // channel: [ versions.yml ]
}
