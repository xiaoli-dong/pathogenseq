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
    ASSEMBYSTATS as STATS_POLYPOLISH;
} from '../../modules/local/assemblystats'

include { FORMATASSEMBLYSTATS as STATS_POLYPOLISH_FORMATASSEMBLYSTATS; } from '../../modules/local/misc'


workflow RUN_POLYPOLISH {

    take:
        illumina_reads
        draft_contigs
    main:
        ch_versions = Channel.empty()

        //draft_contigs.view()
        illumina_reads.map{
            meta, reads ->
                def new_meta = [:]
                new_meta.id = meta.id
                [ new_meta, meta, reads ]

        }.set{
            ch_input_reads
        }

        draft_contigs.map{
            meta, contigs ->
                def new_meta = [:]
                new_meta.id = meta.id
                [ new_meta, meta, contigs ]
        }.set{
            ch_input_contigs
        }

        ch_input_reads.join(ch_input_contigs).multiMap{
            it ->
                reads: [it[1], it[2]]
                draft_contigs: [it[1], it[4]]
        }.set{
            ch_input_all
        }

        BWAMEM2_INDEX(ch_input_all.draft_contigs)

        ch_input_all.reads.join(BWAMEM2_INDEX.out.index).multiMap{
            it ->
                read_1: [it[0], it[1][0]]
                read_2: [it[0], it[1][1]]
                bwa_index: [it[0], it[2]]
        }.set{
            ch_input
        }

        BWAMEM2_MEM_1(ch_input.read_1, ch_input.bwa_index, false)
        BWAMEM2_MEM_2(ch_input.read_2, ch_input.bwa_index, false)

        ch_input_all.draft_contigs.join(BWAMEM2_MEM_1.out.sam).join(BWAMEM2_MEM_2.out.sam).multiMap{
            it ->
                contigs: [it[0], it[1]]
                sam: [it[0], it[2], it[3]]
        }.set{
            ch_input
        }
        POLYPOLISH(ch_input.contigs, ch_input.sam)
        //POLYPOLISH(draft_contigs, sam1.join(sam2))
        ch_versions = ch_versions.mix(POLYPOLISH.out.versions.first())
        contigs = POLYPOLISH.out.contigs
        STATS_POLYPOLISH(contigs)
        STATS_POLYPOLISH_FORMATASSEMBLYSTATS(STATS_POLYPOLISH.out.stats)
        stats = STATS_POLYPOLISH_FORMATASSEMBLYSTATS.out.tsv
        contig_file_ext = ".fasta.gz"

    emit:
        contigs = POLYPOLISH.out.contigs
        contig_file_ext
        versions = ch_versions
        stats

}

include { MASURCA_POLCA } from '../../modules/local/masurca/polca'
include { ASSEMBYSTATS as STATS_POLCA; } from '../../modules/local/assemblystats'
include { FORMATASSEMBLYSTATS as STATS_POLCA_FORMATASSEMBLYSTATS;} from '../../modules/local/misc'

workflow RUN_POLCA{

    take:
        illumina_reads
        draft_contigs
    main:
        ch_versions = Channel.empty()

        //draft_contigs.view()
        illumina_reads.map{
            meta, reads ->
                def new_meta = [:]
                new_meta.id = meta.id
                [ new_meta, meta, reads ]

        }.set{
            ch_input_reads
        }

        draft_contigs.map{
            meta, contigs ->
                def new_meta = [:]
                new_meta.id = meta.id
                [ new_meta, meta, contigs ]
        }.set{
            ch_input_contigs
        }

        ch_input_reads.join(ch_input_contigs).multiMap{
            it ->
                reads: [it[1], it[2]]
                draft_contigs: [it[1], it[4]]
        }.set{
            ch_input
        }

         /* illumina_reads.join(contigs).multiMap{
            it ->
                illumina_reads: [it[0], it[1]]
                contigs: [it[0], it[2]]
        }.set{
            ch_input
        } */

        MASURCA_POLCA(ch_input.reads, ch_input.draft_contigs)
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

include { PYPOLCA} from '../../modules/local/pypolca'
include { ASSEMBYSTATS as STATS_PYPOLCA; } from '../../modules/local/assemblystats'
include {FORMATASSEMBLYSTATS as STATS_PYPOLCA_FORMATASSEMBLYSTATS;} from '../../modules/local/misc'


workflow RUN_PYPOLCA{

    take:
        illumina_reads
        draft_contigs
    main:
        ch_versions = Channel.empty()

        //draft_contigs.view()
        illumina_reads.map{
            meta, reads ->
                def new_meta = [:]
                new_meta.id = meta.id
                [ new_meta, meta, reads ]

        }.set{
            ch_input_reads
        }

        draft_contigs.map{
            meta, contigs ->
                def new_meta = [:]
                new_meta.id = meta.id
                [ new_meta, meta, contigs ]
        }.set{
            ch_input_contigs
        }

        ch_input_reads.join(ch_input_contigs).multiMap{
            it ->
                reads: [it[1], it[2]]
                draft_contigs: [it[1], it[4]]
        }.set{
            ch_input
        }


        PYPOLCA(ch_input.reads, ch_input.draft_contigs)
        contigs = PYPOLCA.out.corrected_contigs
        STATS_PYPOLCA(contigs)
        STATS_PYPOLCA_FORMATASSEMBLYSTATS(STATS_PYPOLCA.out.stats)
        stats = STATS_PYPOLCA_FORMATASSEMBLYSTATS.out.tsv
        ch_versions = ch_versions.mix(PYPOLCA.out.versions.first())
        contig_file_ext = ".fasta.gz"

    emit:
        contigs = PYPOLCA.out.corrected_contigs
        contig_file_ext = ".fasta.gz"
        stats
        versions = ch_versions
}
