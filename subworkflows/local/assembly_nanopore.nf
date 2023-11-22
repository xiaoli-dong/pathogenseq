
include {FLYE} from  '../../modules/nf-core/flye/main'
include {MEDAKA} from  '../../modules/local/medaka/main'

include {
    ASSEMBYSTATS as STATS_MEDAKA;
    ASSEMBYSTATS as STATS_FLYE;

} from '../../modules/local/assemblystats'

include {

    FORMATCSV as STATS_MEDAKA_REFORMAT;
    FORMATCSV as STATS_FLYE_REFORMAT;
    
} from '../../modules/local/formatcsv'


include {RESTARTGENOME } from  '../../modules/local/restartgenome'

workflow ASSEMBLE_NANOPORE {   

    take:
        nanopore_reads
    main:
        ch_versions = Channel.empty()
       
        
        //Flye to be the best-performing bacterial genome assembler in many metrics
        if ( params.nanopore_reads_assembler == 'flye+medaka'){
            
          
            nanopore_reads.multiMap{
                it ->
                nanopore_reads: [it[0], it[1]]
                modeFlag: modeFlag = (
                    it[0].basecaller_mode =~ '_sup_' 
                    || it[0].basecaller_mode =~ '_hac_' 
                ) ? "--nano-hq" : "--nano-raw" 
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
            STATS_FLYE_REFORMAT(STATS_FLYE.out.stats)
            stats = STATS_FLYE_REFORMAT.out.tsv

            //recenter the genome before medaka hopes to fix the termial errors
            RESTARTGENOME(contigs, txt)
            ch_versions = ch_versions.mix(RESTARTGENOME.out.versions.first())

            RESTARTGENOME.out.fasta
                .filter { meta, fasta -> fasta.countFasta() > 0 }
                .set { contigs }
            
            //input = nanopore_reads.join(contigs)
            //nanopore_reads.join(contigs).view()
            nanopore_reads.join(contigs).multiMap{
                it ->
                long_reads_and_assembly: [it[0], it[1], it[2]]
                modeFlag: it[0].basecaller_mode 
            }.set{
                ch_input_medaka
            }
            ch_input_medaka.long_reads_and_assembly.view()
            ch_input_medaka.modeFlag.view()

            MEDAKA(ch_input_medaka.long_reads_and_assembly, ch_input_medaka.modeFlag)

            contigs = MEDAKA.out.assembly
            //SEQKIT_STATS_MEDAKA(contigs)
            STATS_MEDAKA(contigs)
            STATS_MEDAKA_REFORMAT(STATS_MEDAKA.out.stats)
            stats = STATS_MEDAKA_REFORMAT.out.tsv
            contig_file_ext = ".fa.gz"
        } 
        
        
    emit:
        contigs
        contig_file_ext
        versions = ch_versions
        stats
}
