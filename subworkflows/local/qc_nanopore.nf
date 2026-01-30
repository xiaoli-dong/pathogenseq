//nanopore
include {
    NANOPLOT as NANOPLOT_INPUT ;
    NANOPLOT as NANOPLOT_QC
} from '../../modules/nf-core/nanoplot/main'

include { PORECHOP_PORECHOP } from '../../modules/nf-core/porechop/porechop/main.nf'
include { FASTPLONG } from '../../modules/local/fastplong/main.nf'
include { CHOPPER } from '../../modules/local/chopper/main.nf'

include { DEACON_FILTER as DEACON_FILTER_NANOPORE } from '../../modules/local/deacon/filter/main'

include {
    SEQKIT_STATS as SEQKIT_STATS_INPUT_NANOPORE ;
    SEQKIT_STATS as SEQKIT_STATS_PORECHOP ;
    SEQKIT_STATS as SEQKIT_STATS_CHOPPER ;
    SEQKIT_STATS as SEQKIT_STATS_TRIMMED_FASTPLONG ;
    SEQKIT_STATS as SEQKIT_STATS_DEHOST
} from '../../modules/local/seqkit/stats/main'

include {
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_INPUT_NANOPORE ;
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_PORECHOP ;
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_CHOPPER ;
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_DEHOST_NANOPORE
} from '../../modules/nf-core/csvtk/concat/main'

workflow QC_NANOPORE {
    take:
    reads
    adapter_fasta

    main:

    ch_versions = Channel.empty()
    in_format = "tsv"
    out_format = "tsv"

    //reads.view()
    NANOPLOT_INPUT(reads)
    ch_versions = ch_versions.mix(NANOPLOT_INPUT.out.versions.first())
    //reads.view()
    SEQKIT_STATS_INPUT_NANOPORE(reads)
    ch_versions = ch_versions.mix(SEQKIT_STATS_INPUT_NANOPORE.out.versions.first())

    CSVTK_CONCAT_STATS_INPUT_NANOPORE(
        SEQKIT_STATS_INPUT_NANOPORE.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "reads_nanopore.input_seqstats"], files) },
        in_format,
        out_format,
    )

    // QC
    if (params.nanopore_reads_qc_tool == 'fastplong') {
        save_trimmed_fail = false
        save_merged = false
        FASTPLONG(reads, adapter_fasta, save_trimmed_fail)
        ch_versions = ch_versions.mix(FASTPLONG.out.versions.first())
        //get rid of zero size contig file and avoid the downstream crash

        FASTPLONG.out.reads
            .filter { meta, reads -> reads.size() > 0 && reads.countFastq() > 0 }
            .set { trimmed_reads }

        //qc_reads = FASTP.out.reads
        SEQKIT_STATS_TRIMMED_FASTPLONG(trimmed_reads)

        qc_reads = trimmed_reads
        qc_stats = SEQKIT_STATS_TRIMMED_FASTPLONG.out.stats
    }
    else if (params.nanopore_reads_qc_tool == 'porechop_chopper') {
        PORECHOP_PORECHOP(reads)
        ch_versions = ch_versions.mix(PORECHOP_PORECHOP.out.versions.first())
        PORECHOP_PORECHOP.out.reads
            .filter { meta, reads -> reads.size() > 0 && reads.countFastq() > 0 }
            .set { nanopore_reads }

        SEQKIT_STATS_PORECHOP(nanopore_reads)
        porechop_stats = SEQKIT_STATS_PORECHOP.out.stats
        CSVTK_CONCAT_STATS_PORECHOP(
            SEQKIT_STATS_PORECHOP.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "reads_nanopore.porechop_seqstats"], files) },
            in_format,
            out_format,
        )

        CHOPPER(nanopore_reads)
        ch_versions = ch_versions.mix(CHOPPER.out.versions.first())

        CHOPPER.out.fastq
            .filter { meta, reads -> reads.size() > 0 && reads.countFastq() > 0 }
            .set { qc_reads }
        //qc_reads = CHOPPER.out.fastq //gzip compressed
        SEQKIT_STATS_CHOPPER(qc_reads)
        qc_stats = SEQKIT_STATS_CHOPPER.out.stats
    }
    CSVTK_CONCAT_STATS_CHOPPER(
        qc_stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "reads_nanopore.chopper_seqstats"], files) },
        in_format,
        out_format,
    )

    qc_reads.view()
    if (!params.skip_nanopore_dehost) {

        DEACON_FILTER_NANOPORE(
            qc_reads.map { meta, qc_reads -> [meta, params.deacon_refdb, qc_reads] }
        )
        ch_versions = ch_versions.mix(DEACON_FILTER_NANOPORE.out.versions.first())

        qc_reads = DEACON_FILTER_NANOPORE.out.fastq_filtered.filter { meta, reads ->
            def readsList = reads instanceof List ? reads : [reads]
            readsList.size() > 0 && readsList.every { it != null && it.exists() && it.size() > 0 }
        }

        SEQKIT_STATS_DEHOST(qc_reads)
        qc_stats = SEQKIT_STATS_DEHOST.out.stats

        CSVTK_CONCAT_STATS_DEHOST_NANOPORE(
            qc_stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "reads_nanopore.dehost_seqstats"], files) },
            in_format,
            out_format,
        )
    }

    NANOPLOT_QC(qc_reads)

    emit:
    qc_reads
    input_stats = SEQKIT_STATS_INPUT_NANOPORE.out.stats
    qc_stats
    versions = ch_versions
}
