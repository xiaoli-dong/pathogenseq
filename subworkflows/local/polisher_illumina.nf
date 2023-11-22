include {
    MINIMAP2_ALIGN as MINIMAP2_ALIGN1;
    MINIMAP2_ALIGN as MINIMAP2_ALIGN2;
} from '../../modules/local/minimap2/align/main'

include {
    SAMTOOLS_VIEW as SAMTOOLS_VIEW1;
    SAMTOOLS_VIEW as SAMTOOLS_VIEW2;
} from '../../modules/nf-core/samtools/view/main'

include {POLYPOLISH} from '../../modules/local/polypolish'
include {MASURCA_POLCA} from '../../modules/local/masurca/polca'

include{
    SAMTOOLS_SORT as SAMTOOLS_SORT1;
    SAMTOOLS_SORT as SAMTOOLS_SORT2;

} from '../../modules/nf-core/samtools/sort/main'

include {
    ASSEMBYSTATS as STATS_POLYPOLISH;
    ASSEMBYSTATS as STATS_POLCA;

} from '../../modules/local/assemblystats'

include {

    FORMATCSV as STATS_POLYPOLISH_REFORMAT;
    FORMATCSV as STATS_POLCA_REFORMAT;
    
} from '../../modules/local/formatcsv'


workflow RUN_POLYPOLISH {   

    take:
        reads
        contigs
    main:
        ch_versions = Channel.empty()
        
        reads.map { 
            meta, reads -> [ meta, reads[0]] }
            .set { input1 }
        
        reads.map { 
            meta, reads -> [ meta, reads[1]] }
            .set { input2 }
        
        contigs.map{
            meta, contigs -> contigs
        }.set{ fasta }

        bam_format = true
        cigar_paf_format = false
        cigar_bam = false

        MINIMAP2_ALIGN1( input1, fasta, bam_format, cigar_paf_format, cigar_bam)
        MINIMAP2_ALIGN2( input2, fasta, bam_format, cigar_paf_format, cigar_bam)
        SAMTOOLS_SORT1(MINIMAP2_ALIGN1.out.bam)
        SAMTOOLS_SORT2(MINIMAP2_ALIGN2.out.bam)
        bam1 = SAMTOOLS_SORT1.out.bam
        bam2 = SAMTOOLS_SORT2.out.bam
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN1.out.versions.first())

        bam1.map { meta, bam -> [ meta, bam, []] }
            .set { input1 }
        bam2.map { meta, bam -> [ meta, bam, []] }
            .set { input2} 

        SAMTOOLS_VIEW1 ( input1, [[],[]], [] )
        SAMTOOLS_VIEW2 ( input2, [[],[]], [] )
        sam1 = SAMTOOLS_VIEW1.out.sam
        sam2 = SAMTOOLS_VIEW2.out.sam
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW1.out.versions.first())

        POLYPOLISH(contigs, sam1.join(sam2))
        ch_versions = ch_versions.mix(POLYPOLISH.out.versions.first())
        contigs = POLYPOLISH.out.contigs
        STATS_POLYPOLISH(contigs)
        STATS_POLYPOLISH_REFORMAT(STATS_POLYPOLISH.out.stats)
        stats = STATS_POLYPOLISH_REFORMAT.out.tsv

    emit:
        contigs = POLYPOLISH.out.contigs
        versions = ch_versions
        stats

}

workflow RUN_POLCA{   

    take:
        reads
        contigs
    main:
        ch_versions = Channel.empty()
        MASURCA_POLCA(reads, contigs)
        contigs = MASURCA_POLCA.out.contigs
        STATS_POLCA(contigs)
        STATS_POLCA_REFORMAT(STATS_POLCA.out.stats)
        stats = STATS_POLCA_REFORMAT.out.tsv
        ch_versions = ch_versions.mix(MASURCA_POLCA.out.versions.first())
        
    emit:
        contigs = MASURCA_POLCA.out.contigs
        stats
        versions = ch_versions
}