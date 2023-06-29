
include {FLYE} from  '../../modules/nf-core/flye/main'
include {GUNZIP} from  '../../modules/nf-core/gunzip/main'
include {TABIX_BGZIP} from  '../../modules/nf-core/tabix/bgzip/main'
include {MEDAKA} from  '../../modules/nf-core/medaka/main'
include {SEQKIT_STATS as STATS_MEDAKA} from '../../modules/nf-core/seqkit/stats/main'
include {SEQKIT_STATS as STATS_FLYE} from '../../modules/nf-core/seqkit/stats/main'

workflow RUN_ASSEMBLE_LONG {   

    take:
        long_reads
        short_reads
    main:
        ch_versions = Channel.empty()
        //Flye to be the best-performing bacterial genome assembler in many metrics
        if ( params.long_reads_assembler == 'flye+medaka'){
            mode = "--nano-raw"
            FLYE ( long_reads, mode )
            contigs = FLYE.out.fasta
            ch_versions = ch_versions.mix(FLYE.out.versions.first())
            STATS_FLYE(contigs)
            //Medaka cannot accept gzip file as input and it need bgzip files
            GUNZIP(contigs)
            TABIX_BGZIP(GUNZIP.out.gunzip)
            contigs = TABIX_BGZIP.out.output

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