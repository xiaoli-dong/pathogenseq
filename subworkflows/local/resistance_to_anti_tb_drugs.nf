
include { TBPROFILER_CUSTOMIZE as TBPROFILER_PROFILE_ILLUMINA } from '../../modules/local/tbprofiler/customize'    
include { TBPROFILER_CUSTOMIZE as TBPROFILER_PROFILE_NANOPORE } from '../../modules/local/tbprofiler/customize'    

workflow RUN_TBPROFILER_PROFILE_ILLUMINA  {

    take:
        short_reads
        ch_tbprofiler_db

    main:

    ch_versions = Channel.empty()
    
    TBPROFILER_PROFILE_ILLUMINA  (short_reads, ch_tbprofiler_db)
    ch_versions = ch_versions.mix(TBPROFILER_PROFILE_ILLUMINA .out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    bam      = TBPROFILER_PROFILE_ILLUMINA .out.bam         
    csv      = TBPROFILER_PROFILE_ILLUMINA .out.csv       
    json     = TBPROFILER_PROFILE_ILLUMINA .out.json
    txt      =  TBPROFILER_PROFILE_ILLUMINA .out.txt
    vcf      = TBPROFILER_PROFILE_ILLUMINA .out.vcf
    versions = ch_versions                     
}

workflow RUN_TBPROFILER_PROFILE_NANOPORE {

    take:
        long_reads
        ch_tbprofiler_db
    main:

    ch_versions = Channel.empty()
   
    TBPROFILER_PROFILE_NANOPORE (long_reads, ch_tbprofiler_db)
    ch_versions = ch_versions.mix(TBPROFILER_PROFILE_NANOPORE .out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    bam      = TBPROFILER_PROFILE_NANOPORE .out.bam         
    csv      = TBPROFILER_PROFILE_NANOPORE .out.csv       
    json     = TBPROFILER_PROFILE_NANOPORE .out.json
    txt      =  TBPROFILER_PROFILE_NANOPORE .out.txt
    vcf      = TBPROFILER_PROFILE_NANOPORE .out.vcf
    versions = ch_versions                     
}

