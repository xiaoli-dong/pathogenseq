include {
    BWAMEM2_INDEX
} from '../../modules/local/bwamem2/index/main'

include {
    BWAMEM2_MEM as BWAMEM2_MEM_1;
    BWAMEM2_MEM as BWAMEM2_MEM_2;
} from '../../modules/local/bwamem2/mem'


include {
    POLYPOLISH
} from '../../modules/local/polypolish'
include {
    MASURCA_POLCA
} from '../../modules/local/masurca/polca'

include {
    ASSEMBYSTATS as STATS_POLYPOLISH;
    ASSEMBYSTATS as STATS_POLCA;

} from '../../modules/local/assemblystats'

include {

    FORMATASSEMBLYSTATS as STATS_POLYPOLISH_FORMATASSEMBLYSTATS;
    FORMATASSEMBLYSTATS as STATS_POLCA_FORMATASSEMBLYSTATS;
    
} from '../../modules/local/misc'


workflow RUN_POLYPOLISH {   

    take:
        reads
        draft_contigs
    main:
        ch_versions = Channel.empty()
        
        BWAMEM2_INDEX(draft_contigs)

        reads.map { 
            meta, reads -> [ meta, reads[0]] }
            .set { read1 }
        
        reads.map { 
            meta, reads -> [ meta, reads[1]] }
            .set { read2 }
        
        BWAMEM2_MEM_1(read1, BWAMEM2_INDEX.out.index, false)
        BWAMEM2_MEM_2(read2, BWAMEM2_INDEX.out.index, false)

        POLYPOLISH(draft_contigs, BWAMEM2_MEM_1.out.sam.join(BWAMEM2_MEM_2.out.sam))
        //POLYPOLISH(draft_contigs, sam1.join(sam2))
        ch_versions = ch_versions.mix(POLYPOLISH.out.versions.first())
        contigs = POLYPOLISH.out.contigs
        STATS_POLYPOLISH(contigs)
        STATS_POLYPOLISH_FORMATASSEMBLYSTATS(STATS_POLYPOLISH.out.stats)
        stats = STATS_POLYPOLISH_FORMATASSEMBLYSTATS.out.tsv

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
        STATS_POLCA_FORMATASSEMBLYSTATS(STATS_POLCA.out.stats)
        stats = STATS_POLCA_FORMATASSEMBLYSTATS.out.tsv
        ch_versions = ch_versions.mix(MASURCA_POLCA.out.versions.first())
        
    emit:
        contigs = MASURCA_POLCA.out.contigs
        stats
        versions = ch_versions
}