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
} from '../../modules/nf-core/pneumocat/main.nf'
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
     REFORMATEMMTYPERCSV;
} from '../../modules/local/misc'
workflow SPECIAL_TOOLS {   
    take:
        illumina_reads 
        nanopore_reads
        contigs   
    main:
        
        ch_software_versions = Channel.empty()
        // ******* tools for special organisms ************
        //tb-profiler for Mycobacterium tuberculosis
        if(! params.skip_tbprofiler){
            if(illumina_reads){
                TBPROFILER_PROFILE_ILLUMINA(illumina_reads)
                ch_software_versions = ch_software_versions.mix(TBPROFILER_PROFILE_ILLUMINA.out.versions)
                TBPROFILER_COLLATE_ILLUMINA(
                    TBPROFILER_PROFILE_ILLUMINA.out.json.map {meta,json -> json }.collect().map{files -> tuple([id:"tbprofiler"], files)}
                )
            }
            if(nanopore_reads){
                TBPROFILER_PROFILE_NANOPORE(nanopore_reads)
                ch_software_versions = ch_software_versions.mix(TBPROFILER_PROFILE_NANOPORE.out.versions)
                TBPROFILER_COLLATE_NANOPORE(
                    TBPROFILER_PROFILE_NANOPORE.out.json.map {meta,json -> json }.collect().map{files -> tuple([id:"tbprofiler_nanopore"], files)}
                )
            }
        } 
    
        //capsular type to Streptococcus pneumoniae
        if(! params.skip_pneumocat){
            PNEUMOCAT(illumina_reads)
            ch_software_versions = ch_software_versions.mix(PNEUMOCAT.out.versions)
            //PNEUMOCAT.out.xml.map{meta, tsv -> tsv }.collect().map { files -> tuple([id:"pneumocat"], files)}.view()
            COMBINE_XML_PNEUMOCAT(PNEUMOCAT.out.xml.map{meta, tsv -> tsv }.collect().map { files -> tuple([id:"pneumocat"], files)})
            ch_software_versions = ch_software_versions.mix(COMBINE_XML_PNEUMOCAT.out.versions)
        }   
       
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
    
    emit:
        
        versions = ch_software_versions
        
}
