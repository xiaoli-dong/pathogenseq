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
include {
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_ASM;
} from '../../modules/nf-core/csvtk/concat/main'

workflow ASSEMBLE_ILLUMINA {   

    take:
        reads
    main:
        ch_versions = Channel.empty()
        in_format = "tsv"
        out_format = "tsv"
        
        if ( params.illumina_reads_assembler == 'spades' ){
            reads
                .map { meta, reads -> [ meta, reads, [], []] }
                .set { input }
            SPADES (input, [], [])

            //get rid of zero size contig file and avoid the downstream crash
            SPADES.out.contigs
                .filter { meta, contigs -> contigs.countFasta() > 0 }
                .set { contigs }

            //contigs = SPADES.out.contigs
            contig_file_ext = ".fa.gz"
            ch_versions = ch_versions.mix(SPADES.out.versions.first())
            STATS_SPADES(contigs)
            stats = STATS_SPADES.out.stats
            ch_versions = ch_versions.mix(STATS_SPADES.out.versions.first())

        } 
        //default
        else if (params.illumina_reads_assembler == 'skesa' ) {
            SKESA ( reads )
            SKESA.out.contigs
                .filter { meta, contigs -> contigs.countFasta() > 0 }
                .set { contigs }
            contig_file_ext = ".fa.gz"
            //contigs = SKESA.out.contigs
            ch_versions = ch_versions.mix(SKESA.out.versions.first())
            STATS_SKESA(contigs)
            stats = STATS_SKESA.out.stats
            ch_versions = ch_versions.mix(STATS_SKESA.out.versions.first())
        }
        else if (params.illumina_reads_assembler == 'unicycler' ) {
            reads
                .map { meta, reads -> [ meta, reads, [] ] }
                .set { input }
            
            UNICYCLER(input)
            UNICYCLER.out.scaffolds
                .filter { meta, scaffolds -> scaffolds.size() > 0 }
                .set { contigs }

            //contigs = UNICYCLER.out.scaffolds
            contig_file_ext = ".fa.gz"
            ch_versions = ch_versions.mix(UNICYCLER.out.versions.first())
            STATS_UNICYCLER(contigs)
            stats = STATS_UNICYCLER.out.stats
            ch_versions = ch_versions.mix(STATS_UNICYCLER.out.versions.first())
        }
        else if (params.illumina_reads_assembler == 'megahit' ) {
            MEGAHIT ( reads )
            MEGAHIT.out.contigs
                .filter { meta, contigs -> contigs.countFasta() > 0 }
                .set { contigs }
            contig_file_ext = ".fa.gz"
            //contigs = MEGAHIT.out.contigs
            ch_versions = ch_versions.mix(MEGAHIT.out.versions.first())
            STATS_MEGAHIT(contigs)
            stats = STATS_MEGAHIT.out.stats
            ch_versions = ch_versions.mix(STATS_MEGAHIT.out.versions.first())
        }
        CSVTK_CONCAT_STATS_ASM(stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"assembly_stats"], files)}, in_format, out_format ) 
        
    emit:
        contigs
        contig_file_ext
        stats
        versions = ch_versions
        
}
