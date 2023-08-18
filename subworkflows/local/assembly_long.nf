
include {FLYE} from  '../../modules/nf-core/flye/main'
include {GUNZIP} from  '../../modules/nf-core/gunzip/main'
include {TABIX_BGZIP} from  '../../modules/nf-core/tabix/bgzip/main'
include {MEDAKA} from  '../../modules/nf-core/medaka/main'
include {SEQKIT_STATS as STATS_MEDAKA} from '../../modules/nf-core/seqkit/stats/main'
include {SEQKIT_STATS as STATS_FLYE} from '../../modules/nf-core/seqkit/stats/main'
include {RESTARTGENOME } from  '../../modules/local/restartgenome'
//include {DNAAPLER } from  '../../modules/local/dnaapler'

workflow RUN_ASSEMBLE_LONG {   

    take:
        long_reads
        short_reads
    main:
        ch_versions = Channel.empty()
       
        
        //Flye to be the best-performing bacterial genome assembler in many metrics
        if ( params.long_reads_assembler == 'flye+medaka'){
            
            long_reads.multiMap{
                it ->
                long_reads: [it[0], it[1]]
                modeFlag: modeFlag = (it[0].long_mode == 'sup' || it[0].long_mode == 'hac' || it[0].long_mode == 'SUP' || it[0].long_mode == 'HAC') ? "--nano-hq" : "--nano-raw" 
            }.set{
                input
            }
            input.modeFlag.view()
            FLYE(input.long_reads, input.modeFlag)
            FLYE.out.fasta
                .filter { meta, fasta -> fasta.countFasta() > 0 }
                .set { contigs }
            
            FLYE.out.txt
                .filter { meta, txt -> txt.countLines() > 0 }
                .set { txt }

            //contigs = FLYE.out.fasta
            ch_versions = ch_versions.mix(FLYE.out.versions.first())
            STATS_FLYE(contigs)

            //recenter the genome before medaka hopes to fix the termial errors
            RESTARTGENOME(contigs, txt)
            ch_versions = ch_versions.mix(RESTARTGENOME.out.versions.first())

            RESTARTGENOME.out.fasta
                .filter { meta, fasta -> fasta.countFasta() > 0 }
                .set { contigs }
            //contigs = RESTARTGENOME.out.fasta

            //Medaka cannot accept gzip file as input and it need bgzip files
            // GUNZIP(contigs)
            // TABIX_BGZIP(GUNZIP.out.gunzip)
            // contigs = TABIX_BGZIP.out.output

            input = long_reads.join(contigs)
            //file need to be bgzipped
            MEDAKA(input)
            contigs = MEDAKA.out.assembly
            STATS_MEDAKA(contigs)
            stats = STATS_MEDAKA.out.stats
        } 
        
    emit:
        contigs
        versions = ch_versions
        stats
}

/* process seq2fastaFile{
    input:
    tuple val(meta), val(seqid), val(seq)
    output:
    tuple val(meta), val(seqid), path("*.fasta"), emit: fasta

    script:

    """

    echo "${seq}" > ${seqid}.fasta
    

    """
}

workflow run_dnaA_genome{

    take:
        assembly
        
    main:
        ch_versions = Channel.empty()

        assembly.map{
            it ->
            def meta = it[0]
            def arr = it[1].splitFasta( record: [id: true, sequence: true])
            a = []
            arr.each{
                n ->
                a = [meta, n.id, ">"+n.id+"\n" + n.sequence]
            }
            a
        }.set{
            consensus_seq_ch
        }
        consensus_seq_ch.view()
        seq2fastaFile(consensus_seq_ch)
        seq2fastaFile.out.view()
        DNAAPLER(seq2fastaFile.out.fasta)
        fasta = DNAAPLER.out.fasta
        fasta.view()
        

    emit:
        fasta
        versions = ch_versions

        
}
 */