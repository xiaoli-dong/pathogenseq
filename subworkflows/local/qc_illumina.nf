
include {FASTQC as FASTQC_RAW} from     '../../modules/nf-core/fastqc/main'
include {FASTQC as FASTQC_QC} from      '../../modules/nf-core/fastqc/main'
//include { BBMAP_BBDUK} from '../../modules/local/bbmap_bbduk'
include {BBMAP_BBDUK} from      '../../modules/nf-core/bbmap/bbduk/main'
include {FASTP} from '../../modules/nf-core/fastp/main'
include {SEQKIT_STATS as SEQKIT_STATS_RAW} from '../../modules/nf-core/seqkit/stats/main'
include {SEQKIT_STATS as SEQKIT_STATS_QC} from '../../modules/nf-core/seqkit/stats/main'
include {
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_RAW;
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_QC;
} from '../../modules/nf-core/csvtk/concat/main'

workflow QC_ILLUMINA {   

    take:
        reads
    main:
        ch_versions = Channel.empty()
        in_format = "tsv"
        out_format = "tsv"
        qc_reads = reads
        FASTQC_RAW(reads)
        ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())
        SEQKIT_STATS_RAW(reads)
        ch_versions = ch_versions.mix(SEQKIT_STATS_RAW.out.versions.first())
        CSVTK_CONCAT_STATS_RAW(SEQKIT_STATS_RAW.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "illumina_reads_raw.seqstats"], files)}, in_format, out_format )
        
        //default
        if ( params.illumina_reads_qc_tool == 'bbduk' ){
            BBMAP_BBDUK(reads, [])
            ch_versions = ch_versions.mix(BBMAP_BBDUK.out.versions.first())
            //get rid of zero size contig file and avoid the downstream crash
            /* BBMAP_BBDUK.out.reads
                .filter { meta, reads -> reads.countFastq() > 0 }
                .set { qc_reads } */
            qc_reads = BBMAP_BBDUK.out.reads
        }
        else if ( params.illumina_reads_qc_tool == 'fastp'){
            save_trimmed_fail = false
            save_merged       = false
            FASTP ( reads, [], save_trimmed_fail, save_merged )
            ch_versions = ch_versions.mix(FASTP.out.versions.first())
            //get rid of zero size contig file and avoid the downstream crash
            
            FASTP.out.reads
                .filter {meta, reads -> reads[0].size() > 0 && reads[0].countFastq() > 0}
                //.filter { meta, reads -> reads[0].countFastq() > 0 }
                .set { qc_reads }

            //qc_reads = FASTP.out.reads
        }
    
        FASTQC_QC(qc_reads)
       
        SEQKIT_STATS_QC(qc_reads)
        CSVTK_CONCAT_STATS_QC(SEQKIT_STATS_QC.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "illumina_reads_qc.seqstats"], files)}, in_format, out_format )
    
    emit:
        qc_reads
        raw_stats = SEQKIT_STATS_RAW.out.stats
        qc_stats = SEQKIT_STATS_QC.out.stats
        versions = ch_versions
        
}

