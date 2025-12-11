
include {FLYE} from  '../../modules/local/flye/main'


include {
    ASSEMBYSTATS as STATS_MEDAKA;
    ASSEMBYSTATS as STATS_FLYE;

} from '../../modules/local/assemblystats'

include {

    FORMATASSEMBLYSTATS as STATS_MEDAKA_FORMATASSEMBLYSTATS;
    FORMATASSEMBLYSTATS as STATS_FLYE_FORMATASSEMBLYSTATS;

} from '../../modules/local/misc'

include {RESTARTGENOME } from  '../../modules/local/restartgenome'

include { POLISHER_NANOPORE } from './polisher_nanopore'

workflow ASSEMBLE_NANOPORE {

    take:
        nanopore_reads
    main:
        ch_versions = Channel.empty()

        //Flye to be the best-performing bacterial genome assembler in many metrics
        if ( params.nanopore_reads_assembler == 'flye+medaka'){


           /*  nanopore_reads.multiMap{
                it ->
                nanopore_reads: [it[0], it[1]]
                modeFlag: modeFlag = (
                    it[0].basecaller_mode =~ '_sup_'
                    || it[0].basecaller_mode =~ '_hac_'
                ) ? "--nano-hq" : "--nano-raw"
            }.set{
                input
            }
            */
            nanopore_reads.multiMap{
                it ->
                nanopore_reads: [it[0], it[1]]
                modeFlag: modeFlag = (
                    it[0].basecaller_mode =~ 'r941'
                ) ? "--nano-raw" : "--nano-hq"
            }.set{
                input
            }

            //input.modeFlag.view()
            FLYE(input.nanopore_reads, input.modeFlag)
            FLYE.out.fasta
                .filter { meta, fasta -> fasta.countFasta() > 0 }
                .set { contigs }

            FLYE.out.txt
                .filter { meta, txt -> txt.countLines() > 0 }
                .set { txt }

            //contigs = FLYE.out.fasta
            ch_versions = ch_versions.mix(FLYE.out.versions.first())
            //SEQKIT_STATS_FLYE(contigs)
            STATS_FLYE(contigs)
            STATS_FLYE_FORMATASSEMBLYSTATS(STATS_FLYE.out.stats)
            stats = STATS_FLYE_FORMATASSEMBLYSTATS.out.tsv

            //recenter the genome before medaka hopes to fix the termial errors
            if(!params.skip_recenter_genome){
                RESTARTGENOME(contigs, txt)
                ch_versions = ch_versions.mix(RESTARTGENOME.out.versions.first())

                RESTARTGENOME.out.fasta
                    .filter { meta, fasta -> fasta.countFasta() > 0 }
                    .set { contigs }
            }

            POLISHER_NANOPORE(nanopore_reads, contigs)


            contigs = POLISHER_NANOPORE.out.contigs
            STATS = POLISHER_NANOPORE.out.stats
            versions = POLISHER_NANOPORE.out.versions
        }



    emit:
        contigs
        versions = ch_versions
        stats
}
