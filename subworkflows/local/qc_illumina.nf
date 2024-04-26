
include {
    FASTQC as FASTQC_INPUT;
    FASTQC as FASTQC_TRIMMED_BBDUK;
    FASTQC as FASTQC_TRIMMED_FASTP;
    FASTQC as FASTQC_HOSTILE;
} from     '../../modules/nf-core/fastqc/main'

include {
    BBMAP_BBDUK
} from '../../modules/nf-core/bbmap/bbduk/main'
include {
    FASTP
} from '../../modules/local/fastp/main'
include {
    SEQKIT_STATS as SEQKIT_STATS_INPUT_ILLUMINA;
    SEQKIT_STATS as SEQKIT_STATS_TRIMMED_FASTP;
    SEQKIT_STATS as SEQKIT_STATS_TRIMMED_BBDUK;
    SEQKIT_STATS as SEQKIT_STATS_HOSTILE_ILLUMINA;
} from '../../modules/local/seqkit/stats/main'

include {
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_INPUT_ILLUMINA;
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_TRIMMED_ILLUMINA;
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_HOSTILE_ILLUMINA;
} from '../../modules/nf-core/csvtk/concat/main'

include {
    HOSTILE as HOSTILE_ILLUMINA;
} from '../../modules/local/hostile/main'

workflow QC_ILLUMINA {   

    take:
        reads
        adapter_fasta
        hostile_human_ref
    main:
        
        ch_versions = Channel.empty()
        in_format = "tsv"
        out_format = "tsv"
        trimmed_reads = reads
        qc_reads = reads
        
        reads
            .filter {meta, reads -> reads[0].size() > 0 && reads[0].countFastq() > 0}
            .set { reads}

        FASTQC_INPUT(reads)
        ch_versions = ch_versions.mix(FASTQC_INPUT.out.versions.first())
        
        SEQKIT_STATS_INPUT_ILLUMINA(reads)
        ch_versions = ch_versions.mix(SEQKIT_STATS_INPUT_ILLUMINA.out.versions.first())
        CSVTK_CONCAT_STATS_INPUT_ILLUMINA(
            SEQKIT_STATS_INPUT_ILLUMINA.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "reads_illumina.input_seqstats"], files)}, 
            in_format, 
            out_format 
        )
        
        //default
        if ( params.illumina_reads_qc_tool == 'bbduk' ){
            BBMAP_BBDUK(reads, adapter_fasta)
            ch_versions = ch_versions.mix(BBMAP_BBDUK.out.versions.first())
            //get rid of zero size contig file and avoid the downstream crash
            /* BBMAP_BBDUK.out.reads
                .filter { meta, reads -> reads.countFastq() > 0 }
                .set { qc_reads } */
            //qc_reads = BBMAP_BBDUK.out.reads
            trimmed_reads = BBMAP_BBDUK.out.reads
            FASTQC_TRIMMED_BBDUK(trimmed_reads)
            SEQKIT_STATS_TRIMMED_BBDUK(trimmed_reads)

            qc_reads = trimmed_reads
            qc_stats = SEQKIT_STATS_TRIMMED_BBDUK.out.stats
           
        }
        else if ( params.illumina_reads_qc_tool == 'fastp'){
            save_trimmed_fail = false
            save_merged       = false
            FASTP ( reads, adapter_fasta, save_trimmed_fail, save_merged )
            ch_versions = ch_versions.mix(FASTP.out.versions.first())
            //get rid of zero size contig file and avoid the downstream crash
            
            FASTP.out.reads
                .filter {meta, reads -> reads[0].size() > 0 && reads[0].countFastq() > 0}
                .set { trimmed_reads }

            //qc_reads = FASTP.out.reads
            FASTQC_TRIMMED_FASTP(trimmed_reads)
            SEQKIT_STATS_TRIMMED_FASTP(trimmed_reads)

            qc_reads = trimmed_reads
            qc_stats = SEQKIT_STATS_TRIMMED_FASTP.out.stats

        }
        
        CSVTK_CONCAT_STATS_TRIMMED_ILLUMINA(
            qc_stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "reads_illumina.trimmed_seqstats"], files)}, 
            in_format, 
            out_format 
        )
        
       
        //FASTQC_QC(qc_reads)
        if(! params.skip_illumina_dehost){
            HOSTILE_ILLUMINA(trimmed_reads, "bowtie2", hostile_human_ref)
            ch_versions = ch_versions.mix(HOSTILE_ILLUMINA.out.versions.first())

            FASTQC_HOSTILE(HOSTILE_ILLUMINA.out.reads)
            SEQKIT_STATS_HOSTILE_ILLUMINA(HOSTILE_ILLUMINA.out.reads)
            CSVTK_CONCAT_STATS_HOSTILE_ILLUMINA(
                SEQKIT_STATS_HOSTILE_ILLUMINA.out.stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "reads_illumina.dehost_seqstats"], files)}, 
                in_format, 
                out_format 
            )
            qc_reads = HOSTILE_ILLUMINA.out.reads
            qc_stats = SEQKIT_STATS_HOSTILE_ILLUMINA.out.stats
        }

    emit:
        input_stats = SEQKIT_STATS_INPUT_ILLUMINA.out.stats
        qc_reads
        qc_stats
        versions = ch_versions
        
}

