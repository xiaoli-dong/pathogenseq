//https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8380430/
//https://peerj.com/articles/12446.pdf

include { SPADES } from '../../modules/nf-core/spades/main' 
include { UNICYCLER } from '../../modules/nf-core/unicycler/main'
include { MEGAHIT } from '../../modules/nf-core/megahit/main'
include { SKESA } from '../../modules/local/skesa'

include {SEQKIT_STATS as STATS_UNICYCLER} from '../../modules/nf-core/seqkit/stats/main'
include {SEQKIT_STATS as STATS_SKESA} from '../../modules/nf-core/seqkit/stats/main'
include {SEQKIT_STATS as STATS_MEGAHIT} from '../../modules/nf-core/seqkit/stats/main'
include {SEQKIT_STATS as STATS_SPADES} from '../../modules/nf-core/seqkit/stats/main'

include { SHOVILL as SHOVILL_SKESA} from '../../modules/nf-core/shovill/main'
include {SEQKIT_STATS as STATS_SHOVILL_SKESA} from '../../modules/nf-core/seqkit/stats/main'

include { SHOVILL as SHOVILL_MEGAHIT} from '../../modules/nf-core/shovill/main'
include {SEQKIT_STATS as STATS_SHOVILL_MEGAHIT} from '../../modules/nf-core/seqkit/stats/main'

include { SHOVILL as SHOVILL} from '../../modules/nf-core/shovill/main'
include {SEQKIT_STATS as STATS_SHOVILL} from '../../modules/nf-core/seqkit/stats/main'

workflow RUN_ASSEMBLE_SHORT {   

    take:
        reads
    main:
        ch_versions = Channel.empty()

        if ( params.short_reads_assembler == 'spades' ){
            reads
                .map { meta, reads -> [ meta, reads, [], []] }
                .set { input }
            SPADES (input, [], [])
            contigs = SPADES.out.contigs
            ch_versions = ch_versions.mix(SPADES.out.versions.first())
            STATS_SPADES(contigs)
            stats = STATS_SPADES.out.stats
            ch_versions = ch_versions.mix(STATS_SPADES.out.versions.first())

        } 
        //default
        else if (params.short_reads_assembler == 'skesa' ) {
            SKESA ( reads )
            contigs = SKESA.out.contigs
            ch_versions = ch_versions.mix(SKESA.out.versions.first())
            STATS_SKESA(contigs)
            stats = STATS_SKESA.out.stats
            ch_versions = ch_versions.mix(STATS_SKESA.out.versions.first())
        }
        else if (params.short_reads_assembler == 'unicycler' ) {
            reads
                .map { meta, reads -> [ meta, reads, [] ] }
                .set { input }
            
            UNICYCLER(input)
            contigs = UNICYCLER.out.scaffolds
            ch_versions = ch_versions.mix(UNICYCLER.out.versions.first())
            STATS_UNICYCLER(contigs)
            stats = STATS_UNICYCLER.out.stats
            ch_versions = ch_versions.mix(STATS_UNICYCLER.out.versions.first())
        }
        else if (params.short_reads_assembler == 'megahit' ) {
            MEGAHIT ( reads )
            contigs = MEGAHIT.out.contigs
            ch_versions = ch_versions.mix(MEGAHIT.out.versions.first())
            STATS_MEGAHIT(contigs)
            stats = STATS_MEGAHIT.out.stats
            ch_versions = ch_versions.mix(STATS_MEGAHIT.out.versions.first())
        }
        else if (params.short_reads_assembler == 'shovill_skesa' ) {
            SHOVILL_SKESA (reads)
            contigs = SHOVILL_SKESA.out.contigs
            ch_versions = ch_versions.mix(SHOVILL_SKESA.out.versions.first())
            STATS_SHOVILL_SKESA(contigs)
            stats = STATS_SHOVILL_SKESA.out.stats
            ch_versions = ch_versions.mix(STATS_SHOVILL_SKESA.out.versions.first())
        }
        else if (params.short_reads_assembler == 'shovill_megahit' ) {
            SHOVILL_MEGAHIT (reads)
            contigs = SHOVILL_MEGAHIT.out.contigs
            ch_versions = ch_versions.mix(SHOVILL_MEGAHIT.out.versions.first())
            STATS_SHOVILL_MEGAHIT(contigs)
            stats = STATS_SHOVILL_MEGAHIT.out.stats
            ch_versions = ch_versions.mix(STATS_SHOVILL_MEGAHIT.out.versions.first())
        }
        else if (params.short_reads_assembler == 'shovill' ) {
            SHOVILL (reads)
            contigs = SHOVILL.out.contigs
            ch_versions = ch_versions.mix(SHOVILL.out.versions.first())
            STATS_SHOVILL(contigs)
            stats = STATS_SHOVILL.out.stats
            ch_versions = ch_versions.mix(STATS_SHOVILL.out.versions.first())
        }


        
    emit:
        contigs
        stats
        versions = ch_versions
        
}
