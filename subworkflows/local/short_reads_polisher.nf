include {MINIMAP2_ALIGN as MINIMAP2_ALIGN_1} from '../../modules/nf-core/minimap2/align/main'
include {MINIMAP2_ALIGN as MINIMAP2_ALIGN_2} from '../../modules/nf-core/minimap2/align/main'
include {SAMTOOLS_VIEW as SAMTOOLS_VIEW_1} from '../../modules/nf-core/samtools/view/main'
include {SAMTOOLS_VIEW as SAMTOOLS_VIEW_2} from '../../modules/nf-core/samtools/view/main'
include {POLYPOLISH} from '../../modules/local/polypolish'
include {MASURCA_POLCA} from '../../modules/local/masurca/polca'
include {SEQKIT_STATS as POLYPOLISH_STATS} from '../../modules/nf-core/seqkit/stats/main'
include {SEQKIT_STATS as POLCA_STATS} from '../../modules/nf-core/seqkit/stats/main'


workflow RUN_POLYPOLISH {   

    take:
        reads
        assembly
    main:
        ch_versions = Channel.empty()
        
        reads.map { 
            meta, reads -> [ meta, reads[0]] }
            .set { input1 }
        
        reads.map { 
            meta, reads -> [ meta, reads[1]] }
            .set { input2 }
        
        assembly.map{
            meta, assembly -> assembly
        }.set{ fasta }

        bam_format = true
        cigar_paf_format = false
        cigar_bam = false

        MINIMAP2_ALIGN_1( input1, fasta, bam_format, cigar_paf_format, cigar_bam)
        MINIMAP2_ALIGN_2( input2, fasta, bam_format, cigar_paf_format, cigar_bam)
        bam_1 = MINIMAP2_ALIGN_1.out.bam
        bam_2 = MINIMAP2_ALIGN_2.out.bam
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN_1.out.versions.first())

        bam_1.map { meta, bam -> [ meta, bam, []] }
            .set { input1 }
        bam_2.map { meta, bam -> [ meta, bam, []] }
            .set { input2} 

        SAMTOOLS_VIEW_1 ( input1, [[],[]], [] )
        SAMTOOLS_VIEW_2 ( input2, [[],[]], [] )
        sam_1 = SAMTOOLS_VIEW_1.out.sam
        sam_2 = SAMTOOLS_VIEW_2.out.sam
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW_1.out.versions.first())

        POLYPOLISH(assembly, sam_1.join(sam_2))
        ch_versions = ch_versions.mix(POLYPOLISH.out.versions.first())
        assembly = POLYPOLISH.out.assembly 
        POLYPOLISH_STATS(assembly)

    emit:
        assembly = POLYPOLISH.out.assembly
        versions = ch_versions
        stats = POLYPOLISH_STATS.out.stats

}

workflow RUN_POLCA{   

    take:
        reads
        assembly
    main:
        ch_versions = Channel.empty()
        MASURCA_POLCA(reads, assembly)
        contigs = MASURCA_POLCA.out.assembly
        //stats = POLCA_STATS(contigs)

        ch_versions = ch_versions.mix(MASURCA_POLCA.out.versions.first())
        
    emit:
        contigs
        //stats
        versions = ch_versions
}


