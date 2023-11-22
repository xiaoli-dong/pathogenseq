//nanopore
include {
    NANOPLOT as NANOPLOT_INPUT;
    NANOPLOT as NANOPLOT_QC;
} from '../../modules/nf-core/nanoplot/main' 

include {PORECHOP_PORECHOP}  from '../../modules/nf-core/porechop/porechop/main.nf'  
include {CHOPPER}  from '../../modules/local/chopper/main.nf'  

include {
    SEQKIT_STATS as SEQKIT_STATS_INPUT
    SEQKIT_STATS as SEQKIT_STATS_PORECHOP;
    SEQKIT_STATS as SEQKIT_STATS_CHOPPER;
} from '../../modules/local/seqkit/stats/main'

include {
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_INPUT;
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_PORECHOP;
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_CHOPPER;
} from '../../modules/nf-core/csvtk/concat/main'

workflow QC_NANOPORE {

    take:
        reads
    main:

        ch_versions = Channel.empty()
        //reads.view()
        NANOPLOT_INPUT(reads)
        ch_versions = ch_versions.mix(NANOPLOT_INPUT.out.versions.first())
        //reads.view()
        SEQKIT_STATS_INPUT(reads)
        ch_versions = ch_versions.mix(SEQKIT_STATS_INPUT.out.versions.first())
        
        // QC
        PORECHOP_PORECHOP(reads)
        SEQKIT_STATS_PORECHOP(PORECHOP_PORECHOP.out.reads)

        ch_versions = ch_versions.mix(PORECHOP_PORECHOP.out.versions.first())
        CHOPPER(PORECHOP_PORECHOP.out.reads)
        SEQKIT_STATS_CHOPPER(CHOPPER.out.fastq)

        qc_reads = CHOPPER.out.fastq //gzip compressed
        ch_versions = ch_versions.mix(CHOPPER.out.versions.first())
        
        
        NANOPLOT_QC(qc_reads)
        
        in_format = "tsv"
        out_format = "tsv"
        CSVTK_CONCAT_STATS_INPUT(SEQKIT_STATS_INPUT.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"nanopore_reads.input_seqstats"], files)}, in_format, out_format )
        CSVTK_CONCAT_STATS_PORECHOP(SEQKIT_STATS_PORECHOP.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"nanopore_reads.porechop_seqstats"], files)}, in_format, out_format ) 
        CSVTK_CONCAT_STATS_CHOPPER(SEQKIT_STATS_CHOPPER.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"nanopore_reads.chopper_seqstats"], files)}, in_format, out_format ) 
    emit:
        qc_reads
        input_stats = SEQKIT_STATS_INPUT.out.stats
        porechop_stats = SEQKIT_STATS_PORECHOP.out.stats
        qc_stats = SEQKIT_STATS_CHOPPER.out.stats
        versions = ch_versions

}
