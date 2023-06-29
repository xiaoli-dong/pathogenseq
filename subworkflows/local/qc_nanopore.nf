//nanopore
include {NANOPLOT as NANOPLOT_RAW} from '../../modules/nf-core/nanoplot/main' 
include {NANOPLOT as NANOPLOT_QC}  from '../../modules/nf-core/nanoplot/main'
include {PORECHOP_PORECHOP}  from '../../modules/nf-core/porechop/porechop/main.nf'  
include {SEQKIT_STATS as SEQKIT_STATS_LONG_RAW} from '../../modules/nf-core/seqkit/stats/main'
include {SEQKIT_STATS as SEQKIT_STATS_LONG_QC} from '../../modules/nf-core/seqkit/stats/main'

workflow RUN_NANOPORE_QC {

    take:
        reads
    main:

        ch_versions = Channel.empty()
        //reads.view()
        NANOPLOT_RAW(reads)
        ch_versions = ch_versions.mix(NANOPLOT_RAW.out.versions.first())
        //reads.view()
        SEQKIT_STATS_LONG_RAW(reads)
        ch_versions = ch_versions.mix(SEQKIT_STATS_LONG_RAW.out.versions.first())
        
        // QC
        PORECHOP_PORECHOP(reads)
        qc_reads = PORECHOP_PORECHOP.out.reads
        ch_versions = ch_versions.mix(PORECHOP_PORECHOP.out.versions.first())
        SEQKIT_STATS_LONG_QC(qc_reads)

        //qc_reads.view()
        NANOPLOT_QC(qc_reads)


    emit:
        qc_reads
        raw_stats = SEQKIT_STATS_LONG_RAW.out.stats
        qc_stats = SEQKIT_STATS_LONG_QC.out.stats
        versions = ch_versions

}
