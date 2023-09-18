
include {FLYE} from  '../../modules/nf-core/flye/main'
include {MEDAKA} from  '../../modules/nf-core/medaka/main'
include {SEQKIT_STATS as SEQKIT_STATS_MEDAKA} from '../../modules/nf-core/seqkit/stats/main'
include {SEQKIT_STATS as SEQKIT_STATS_FLYE} from '../../modules/nf-core/seqkit/stats/main'
include {RESTARTGENOME } from  '../../modules/local/restartgenome'
//include {DNAAPLER } from  '../../modules/local/dnaapler'

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
                modeFlag: modeFlag = (it[0].basecaller_mode == 'sup' || it[0].basecaller_mode == 'hac' || it[0].basecaller_mode == 'SUP' || it[0].basecaller_mode == 'HAC') ? "--nano-hq" : "--nano-raw" 
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
            SEQKIT_STATS_FLYE(contigs)

            //recenter the genome before medaka hopes to fix the termial errors
            RESTARTGENOME(contigs, txt)
            ch_versions = ch_versions.mix(RESTARTGENOME.out.versions.first())

            RESTARTGENOME.out.fasta
                .filter { meta, fasta -> fasta.countFasta() > 0 }
                .set { contigs }
            
            input = nanopore_reads.join(contigs)
            MEDAKA(input)
            contigs = MEDAKA.out.assembly
            SEQKIT_STATS_MEDAKA(contigs)
            stats = SEQKIT_STATS_MEDAKA.out.stats
            contig_file_ext = ".fa.gz"
        } 
        
        
    emit:
        contigs
        contig_file_ext
        versions = ch_versions
        stats
}
