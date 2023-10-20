include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_DEPTH_NANOPORE } from '../../modules/local/minimap2/align/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_DEPTH_NANOPORE } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEPTH_NANOPORE } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_DEPTH_NANOPORE } from '../../modules/nf-core/samtools/coverage/main'
include { MAPPINGREPORT as MAPPINGREPORT_NANOPORE } from '../../modules/local/mappingreport'
workflow DEPTH_NANOPORE {   

    take:
        reads
        contigs
    main:
        ch_versions = Channel.empty()
        
        contigs.map{
            meta, contigs -> contigs
        }.set{ fasta }

        bam_format = true
        cigar_paf_format = false
        cigar_bam = false

        MINIMAP2_ALIGN_DEPTH_NANOPORE( reads, fasta, bam_format, cigar_paf_format, cigar_bam)
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN_DEPTH_NANOPORE.out.versions.first())

        SAMTOOLS_SORT_DEPTH_NANOPORE(MINIMAP2_ALIGN_DEPTH_NANOPORE.out.bam)
        SAMTOOLS_INDEX_DEPTH_NANOPORE(SAMTOOLS_SORT_DEPTH_NANOPORE.out.bam)
        SAMTOOLS_COVERAGE_DEPTH_NANOPORE(SAMTOOLS_SORT_DEPTH_NANOPORE.out.bam.join(SAMTOOLS_INDEX_DEPTH_NANOPORE.out.bai))
        ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE_DEPTH_NANOPORE.out.versions.first())
        MAPPINGREPORT_NANOPORE(SAMTOOLS_COVERAGE_DEPTH_NANOPORE.out.coverage)

    emit:

        contig_coverage = SAMTOOLS_COVERAGE_DEPTH_NANOPORE.out.coverage
        sample_coverage = MAPPINGREPORT_NANOPORE.out.tsv
        versions = ch_versions
       

}
