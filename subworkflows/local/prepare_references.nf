//
// prepare databases
//

include { AMRFINDERPLUS_UPDATE          } from '../../modules/nf-core/amrfinderplus/update/main'

workflow PREPARE_REFERENCES {
    main:

    ch_versions = Channel.empty()

     //
    // hostile reference
    //
    ch_hostile_ref_minimap2 = Channel.empty()
    if (params.hostile_human_ref_minimap2) {
        ch_hostile_ref_minimap2 = Channel.value(file(params.hostile_human_ref_minimap2))
    }
    else {
        log.error "Please specify a valid hostime human reference file for minimap2"
        System.exit(1)
    }

    ch_hostile_ref_bowtie2 = Channel.empty()
    if (params.hostile_human_ref_bowtie2) {
        ch_hostile_ref_bowtie2 = Channel.value(file(params.hostile_human_ref_bowtie2))
    }
    else {
        log.error "Please specify a valid hostime human reference file from bowtie2"
        System.exit(1)
    }

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

    //
    // Prepare amrfinderplus database
    //
    ch_gambit_db = Channel.empty()
    if (!params.skip_gambit) {
        if (params.gambit_db) {
            ch_gambit_db = Channel.value(file(params.gambit_db))
        }
        else {
            log.error "Please specify a valid  bakta database."
            System.exit(1)
        }
    }

    //
    // Prepare GBS-SBG database
    //
    ch_gbssbg_db = Channel.empty()
    if (!params.skip_gbssbg) {
        if (params.gbssbg_db) {
            ch_gbssbg_db = Channel.value(file(params.gbssbg_db))
        }
        else {
            log.error "Please specify a valid  gbs-sbg database."
            System.exit(1)
        }
    }

    emit:
        ch_hostile_ref_bowtie2
        ch_hostile_ref_minimap2 
        ch_kraken2_db                    // path: kraken2_db/
        ch_checkm2_db   //path to checkm2 db
        ch_bakta_db     // path to bakta database 
        ch_amrfinderplus_db // path: path to amrfinderplus_db
        ch_gambit_db // path: path to gambit db
        ch_gbssbg_db // path: path to gbs-sbg database
        versions             = ch_versions             // channel: [ versions.yml ]
}
