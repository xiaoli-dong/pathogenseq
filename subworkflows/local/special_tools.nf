include {
    TBPROFILER_PROFILE as TBPROFILER_PROFILE_ILLUMINA;
    TBPROFILER_PROFILE as TBPROFILER_PROFILE_NANOPORE
} from '../../modules/local/tbprofiler/profile'
include {
    TBPROFILER_COLLATE as TBPROFILER_COLLATE_ILLUMINA;
     TBPROFILER_COLLATE as TBPROFILER_COLLATE_NANOPORE;
} from '../../modules/local/tbprofiler/collate'

include {
    PNEUMOCAT
} from '../../modules/local/pneumocat/main.nf'
include {
    COMBINE_XML as COMBINE_XML_PNEUMOCAT;
} from '../../modules/local/combine/xml'
include {
    EMMTYPER
} from '../../modules/local/emmtyper/main.nf'
include {
     CSVTK_CONCAT as CSVTK_CONCAT_EMMTYPER;
    
} from '../../modules/nf-core/csvtk/concat/main'

include {
     CSVTK_CONCAT as CSVTK_CONCAT_GBSSBG;
} from '../../modules/local/csvtk/concat/main'

include {
     REFORMATEMMTYPERCSV;
} from '../../modules/local/misc'

include {
    GBS_SBG;
} from '../../modules/local/gbs/sbg'

include {
    PREPARE_REFERENCES
} from '../../subworkflows/local/prepare_references'

workflow SPECIAL_TOOLS_BASED_ON_ILLUMINA {   
    take:
        illumina_reads 
    main:
        
        ch_software_versions = Channel.empty()
        // ******* tools for special organisms ************
        //tb-profiler for Mycobacterium tuberculosis
        if(! params.skip_tbprofiler){
            //if(!illumina_reads.ifEmpty(null)){
                TBPROFILER_PROFILE_ILLUMINA(illumina_reads)
                ch_software_versions = ch_software_versions.mix(TBPROFILER_PROFILE_ILLUMINA.out.versions)
                TBPROFILER_COLLATE_ILLUMINA(
                    TBPROFILER_PROFILE_ILLUMINA.out.json.map {meta,json -> json }.collect().map{files -> tuple([id:"tbprofiler"], files)}
                )
            //}
        } 

        //capsular type to Streptococcus pneumoniae
        //if(! params.skip_pneumocat &&  !illumina_reads.ifEmpty(null)){
        if(! params.skip_pneumocat){
            PNEUMOCAT(illumina_reads)
            ch_software_versions = ch_software_versions.mix(PNEUMOCAT.out.versions)
            COMBINE_XML_PNEUMOCAT(PNEUMOCAT.out.results.map{meta, tsv -> tsv }.collect().map { files -> tuple([id:"pneumocat"], files)})
            ch_software_versions = ch_software_versions.mix(COMBINE_XML_PNEUMOCAT.out.versions)
        }   
       
    
    emit:
        versions = ch_software_versions
        
}

workflow SPECIAL_TOOLS_BASED_ON_NANOPORE {   
    take:
        nanopore_reads
    main:
        
        ch_software_versions = Channel.empty()
        // ******* tools for special organisms ************
        //tb-profiler for Mycobacterium tuberculosis
        if(! params.skip_tbprofiler){
            
            //if(!nanopore_reads.ifEmpty(null)){
                TBPROFILER_PROFILE_NANOPORE(nanopore_reads)
                ch_software_versions = ch_software_versions.mix(TBPROFILER_PROFILE_NANOPORE.out.versions)
                TBPROFILER_COLLATE_NANOPORE(
                    TBPROFILER_PROFILE_NANOPORE.out.json.map {meta,json -> json }.collect().map{files -> tuple([id:"tbprofiler_nanopore"], files)}
                )
            }
        //} 
    
    emit:
        versions = ch_software_versions
        
}

workflow SPECIAL_TOOLS_BASED_ON_CONTIGS {   
    take:
        contigs   
    main:
        
        ch_software_versions = Channel.empty()
       
        PREPARE_REFERENCES ()
        ch_software_versions = ch_software_versions.mix(PREPARE_REFERENCES.out.versions)
       
        //emm-typing of Streptococcus pyogenes
        if(! params.skip_emmtyper){
            EMMTYPER(contigs)
            ch_software_versions = ch_software_versions.mix(EMMTYPER.out.versions)
            REFORMATEMMTYPERCSV(EMMTYPER.out.tsv)
            ch_software_versions = ch_software_versions.mix(REFORMATEMMTYPERCSV.out.versions)
            CSVTK_CONCAT_EMMTYPER(
                REFORMATEMMTYPERCSV.out.csv.map {meta, csv -> csv }.collect().map { files -> tuple([id:"emmtyper"], files)}, 
                "csv", 
                "csv" 
            )
            ch_software_versions = ch_software_versions.mix(CSVTK_CONCAT_EMMTYPER.out.versions)
        }  

         // serotyping GBS: Streptococcus agalactiae
        if(! params.skip_gbssbg){
            GBS_SBG(contigs, PREPARE_REFERENCES.out.ch_gbssbg_db)
            ch_software_versions = ch_software_versions.mix(GBS_SBG.out.versions)
            GBS_SBG.out.tsv.map {meta, tsv -> tsv }.collect()//.view()
            CSVTK_CONCAT_GBSSBG(
                GBS_SBG.out.tsv.map {meta, tsv -> tsv }.collect().map { files -> tuple([id:"gbs-sbg"], files)}, 
                "tsv", 
                "tsv" 
            )
            ch_software_versions = ch_software_versions.mix(CSVTK_CONCAT_GBSSBG.out.versions)
        }  
    
    emit:
        //emm-typing = REFORMATEMMTYPERCSV.out.csv //csv format emmtyper verbose format
        versions = ch_software_versions
        
}

