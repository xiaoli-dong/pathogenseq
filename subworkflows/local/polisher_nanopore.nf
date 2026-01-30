include { MEDAKA_CONSENSUS } from '../../modules/local/medaka/consensus'
include { ASSEMBYSTATS as STATS_MEDAKA } from '../../modules/local/assemblystats'
include { FORMATASSEMBLYSTATS as STATS_MEDAKA_FORMATASSEMBLYSTATS } from '../../modules/local/misc'
workflow POLISHER_NANOPORE {
    take:
    nanopore_reads
    contigs

    main:
    ch_versions = channel.empty()
    stats = channel.empty()

    /* nanopore_reads
        .toList()
        .then {
            def needs_polish = items.findAll { it[0].basecaller_mode != 'NA' }
            def skip_polish = items.findAll { it[0].basecaller_mode == 'NA' }

            if (needs_polish) {
                if (params.nanopore_reads_polisher == 'medaka') {

                    nanopore_reads
                        .join(contigs)
                        .multiMap { meta, reads, contigs ->
                            long_reads_and_assembly: [meta, reads, contigs]
                            modeFlag: meta.basecaller_mode
                        }
                        .set {
                            ch_input_medaka
                        }
                    MEDAKA_CONSENSUS(ch_input_medaka.long_reads_and_assembly, ch_input_medaka.modeFlag)
                    ch_versions = ch_versions.mix(MEDAKA_CONSENSUS.out.versions.first())

                    contigs = MEDAKA_CONSENSUS.out.assembly
                    STATS_MEDAKA(contigs)
                    STATS_MEDAKA_FORMATASSEMBLYSTATS(STATS_MEDAKA.out.stats)

                    stats = STATS_MEDAKA_FORMATASSEMBLYSTATS.out.tsv
                }
            }
        } */

    /*
    * 1. Split nanopore reads by whether polishing is needed
    */
    nanopore_reads
        .branch {
            needs_polish: it[0].basecaller_mode != 'NA'
            skip_polish:  it[0].basecaller_mode == 'NA'
        }
        .set { ch_nanopore_reads }

    /*
    * 2. Prepare Medaka input ONLY for samples that need polishing
    */
    ch_nanopore_reads.needs_polish
        .join(contigs)
        .multiMap { meta, reads, contigs_fa ->
            long_reads_and_assembly: [meta, reads, contigs_fa]
            modeFlag: meta.basecaller_mode
        }
        .set { ch_input_medaka }

    /*
    * 3. Gate Medaka execution (no samples → empty channel)
    */
    ch_input_medaka = (
        params.nanopore_reads_polisher == 'medaka'
            ? ch_input_medaka
            : Channel.empty()
    )

    /*
    * 4. Run Medaka (process is ALWAYS declared)
    */
    MEDAKA_CONSENSUS(
        ch_input_medaka.long_reads_and_assembly,
        ch_input_medaka.modeFlag
    )

    /*
    * 5. Update versions
    */
    ch_versions = ch_versions.mix(MEDAKA_CONSENSUS.out.versions)

    /*
    * 6. Merge polished + unpolished contigs
    *    (polished overrides originals for matching samples)
    */
    contigs = contigs.mix(MEDAKA_CONSENSUS.out.assembly)

    /*
    * 7. Stats on final contigs
    */
    STATS_MEDAKA(contigs)
    STATS_MEDAKA_FORMATASSEMBLYSTATS(STATS_MEDAKA.out.stats)

    stats = STATS_MEDAKA_FORMATASSEMBLYSTATS.out.tsv



    emit:
    contigs
    versions = ch_versions
    stats
}
