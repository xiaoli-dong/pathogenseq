
include {FASTQC as FASTQC_RAW} from     '../../modules/nf-core/fastqc/main'
include {FASTQC as FASTQC_QC} from      '../../modules/nf-core/fastqc/main'
//include { BBMAP_BBDUK} from '../../modules/local/bbmap_bbduk'
include {BBMAP_BBDUK} from      '../../modules/nf-core/bbmap/bbduk/main'
include {FASTP} from '../../modules/nf-core/fastp/main'
include {SEQKIT_STATS as SEQKIT_STATS_SHORT_RAW} from '../../modules/nf-core/seqkit/stats/main'
include {SEQKIT_STATS as SEQKIT_STATS_SHORT_QC} from '../../modules/nf-core/seqkit/stats/main'

workflow RUN_ILLUMINA_QC {   

    take:
        reads
    main:
        ch_versions = Channel.empty()
        qc_reads = reads
        FASTQC_RAW(reads)
        ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())
        SEQKIT_STATS_SHORT_RAW(reads)
        ch_versions = ch_versions.mix(SEQKIT_STATS_SHORT_RAW.out.versions.first())
        
        //default
        if ( params.short_reads_qc_tool == 'bbduk' ){
            BBMAP_BBDUK(reads, [])
            ch_versions = ch_versions.mix(BBMAP_BBDUK.out.versions.first())
            qc_reads = BBMAP_BBDUK.out.reads
        }
        else if ( params.short_reads_qc_tool == 'fastp'){
            save_trimmed_fail = false
            save_merged       = false
            FASTP ( reads, [], save_trimmed_fail, save_merged )
            ch_versions = ch_versions.mix(FASTP.out.versions.first())
            qc_reads = FASTP.out.reads
        }

        FASTQC_QC(qc_reads)
       
        SEQKIT_STATS_SHORT_QC(qc_reads)
        
    emit:
        qc_reads
        raw_stats = SEQKIT_STATS_SHORT_RAW.out.stats
        qc_stats = SEQKIT_STATS_SHORT_QC.out.stats
        versions = ch_versions
        
}

