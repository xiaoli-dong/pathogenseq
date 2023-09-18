//nanopore
include {NANOPLOT as NANOPLOT_RAW} from '../../modules/nf-core/nanoplot/main' 
include {NANOPLOT as NANOPLOT_QC}  from '../../modules/nf-core/nanoplot/main'
include {PORECHOP_PORECHOP}  from '../../modules/nf-core/porechop/porechop/main.nf'  
include {SEQKIT_STATS as SEQKIT_STATS_RAW} from '../../modules/nf-core/seqkit/stats/main'
include {SEQKIT_STATS as SEQKIT_STATS_QC} from '../../modules/nf-core/seqkit/stats/main'
include {
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_RAW;
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_QC;
} from '../../modules/nf-core/csvtk/concat/main'

workflow QC_NANOPORE {

    take:
        reads
    main:

        ch_versions = Channel.empty()
        //reads.view()
        NANOPLOT_RAW(reads)
        ch_versions = ch_versions.mix(NANOPLOT_RAW.out.versions.first())
        //reads.view()
        SEQKIT_STATS_RAW(reads)
        ch_versions = ch_versions.mix(SEQKIT_STATS_RAW.out.versions.first())
        
        // QC
        PORECHOP_PORECHOP(reads)
        qc_reads = PORECHOP_PORECHOP.out.reads
        ch_versions = ch_versions.mix(PORECHOP_PORECHOP.out.versions.first())
        SEQKIT_STATS_QC(qc_reads)

        //qc_reads.view()
        NANOPLOT_QC(qc_reads)
        
        in_format = "tsv"
        out_format = "tsv"
        CSVTK_CONCAT_STATS_RAW(SEQKIT_STATS_RAW.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"nanopore_reads_raw.seqstats"], files)}, in_format, out_format )
        CSVTK_CONCAT_STATS_QC(SEQKIT_STATS_QC.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"nanopore_reads_qc.seqstats"], files)}, in_format, out_format ) 

    emit:
        qc_reads
        raw_stats = SEQKIT_STATS_RAW.out.stats
        qc_stats = SEQKIT_STATS_QC.out.stats
        versions = ch_versions

}
