include {MEDAKA_CONSENSUS} from  '../../modules/local/medaka/consensus'
include { ASSEMBYSTATS as STATS_MEDAKA } from '../../modules/local/assemblystats'
include {FORMATASSEMBLYSTATS as STATS_MEDAKA_FORMATASSEMBLYSTATS} from '../../modules/local/misc'
workflow POLISHER_NANOPORE {
    take:
        nanopore_reads
        contigs
        contig_file_ext
    main:
        ch_versions = Channel.empty()


        
        stats = Channel.empty()
        nanopore_reads.view()
        nanopore_reads.toList().then {
            def needs_polish = items.findAll { it[0].basecaller_mode != 'NA' }
            def skip_polish  = items.findAll { it[0].basecaller_mode == 'NA' }
        
            if(needs_polish){
                if(params.nanopore_reads_polisher == 'medaka'){
                    
                    nanopore_reads.join(contigs).multiMap{
                        meta, reads, contigs  ->
                        long_reads_and_assembly: [meta, reads, contigs]
                        modeFlag: meta.basecaller_mode
                    }.set{
                        ch_input_medaka
                    }
                    MEDAKA_CONSENSUS(ch_input_medaka.long_reads_and_assembly, ch_input_medaka.modeFlag)
                    ch_versions = ch_versions.mix(MEDAKA_CONSENSUS.out.versions.first())

                    contigs = MEDAKA_CONSENSUS.out.assembly
                    STATS_MEDAKA(contigs)
                    STATS_MEDAKA_FORMATASSEMBLYSTATS(STATS_MEDAKA.out.stats)

                    stats = STATS_MEDAKA_FORMATASSEMBLYSTATS.out.tsv
                    contig_file_ext = ".fasta.gz"
                }
                else if(params.nanopore_reads_polisher == 'dorado'){

                }
            }
        }
    emit:
        contigs
        contig_file_ext
        versions = ch_versions
        stats

}
