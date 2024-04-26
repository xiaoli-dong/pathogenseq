/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def valid_params = [
    nanopore_reads_assembler : ['flye+medaka']
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowNanopore.initialise(params, log,  valid_params)
// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
//def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
def checkPathParamList = [ params.input, params.multiqc_config]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPARE_REFERENCES            } from '../subworkflows/local/prepare_references'
include { INPUT_CHECK                   } from '../subworkflows/local/input_check'
include { QC_ILLUMINA                   } from '../subworkflows/local/qc_illumina'
include { QC_NANOPORE                   } from '../subworkflows/local/qc_nanopore'
include { ASSEMBLE_ILLUMINA             } from '../subworkflows/local/assembly_illumina'
include { ANNOTATION                    } from '../subworkflows/local/annotation'
include { ASSEMBLE_NANOPORE             } from '../subworkflows/local/assembly_nanopore'

include { RUN_POLYPOLISH; RUN_POLCA;    } from '../subworkflows/local/polisher_illumina'
include { DEPTH_ILLUMINA    }       from '../subworkflows/local/depth_illumina'
include { DEPTH_NANOPORE   }       from '../subworkflows/local/depth_nanopore'
include { 
   
    SPECIAL_TOOLS_BASED_ON_CONTIGS;
    SPECIAL_TOOLS_BASED_ON_ILLUMINA;
    SPECIAL_TOOLS_BASED_ON_NANOPORE;
} from '../subworkflows/local/special_tools'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { 
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_ASM; 
    CSVTK_CONCAT as CSVTK_CONCAT_DEPTH_NANOPORE;
    CSVTK_CONCAT as CSVTK_CONCAT_DEPTH_ILLUMINA;  
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_NOT_ASSEMBLED;
} from '../modules/nf-core/csvtk/concat/main'

include { 
    KRAKEN2_KRAKEN2 as KRAKEN2_KRAKEN2_ILLUMINA;
    KRAKEN2_KRAKEN2 as KRAKEN2_KRAKEN2_NANOPORE;
} from '../modules/nf-core/kraken2/kraken2/main' 
include { 
    KRAKENTOOLS_COMBINEKREPORTS as KRAKENTOOLS_COMBINEKREPORTS_ILLUMINA;
    KRAKENTOOLS_COMBINEKREPORTS as KRAKENTOOLS_COMBINEKREPORTS_NANOPORE;
} from '../modules/nf-core/krakentools/combinekreports/main.nf'

//MODULES: local modules
include {CHECKM2_PREDICT} from '../modules/local/checkm2/predict.nf'
include { 
    GAMBIT_QUERY as GAMBIT_QUERY_COLLECT;
    GAMBIT_QUERY as GAMBIT_QUERY;
 }    from '../modules/local/gambit/query/main'
include { GAMBIT_TREE       }       from '../modules/local/gambit/tree/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow NANOPORE {

    ch_software_versions = Channel.empty()

    //
    // SUBWORKFLOW: prepare reference databases ...
    //
    PREPARE_REFERENCES ()
    ch_software_versions = ch_software_versions.mix(PREPARE_REFERENCES.out.versions)

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_software_versions = ch_software_versions.mix(INPUT_CHECK.out.versions)

    reads = INPUT_CHECK.out.reads
    illumina_reads = INPUT_CHECK.out.shortreads
    nanopore_reads = INPUT_CHECK.out.longreads
    in_format = "tsv"
    out_format = "tsv"
    contig_file_ext = ".fa.gz"

   
    
    if(!params.skip_illumina_reads_qc){
      
        QC_ILLUMINA(
            illumina_reads.filter {meta, reads -> reads[0].size() > 0 && reads[0].countFastq() > 0},
            [],
            PREPARE_REFERENCES.out.ch_hostile_ref_bowtie2
        )
        ch_software_versions = ch_software_versions.mix(QC_ILLUMINA.out.versions)
        
        //get rid of zero size contig file and avoid the downstream crash
        QC_ILLUMINA.out.qc_reads
            .filter {meta, reads -> reads[0].size() > 0 && reads[0].countFastq() > 0}
            .set { illumina_reads }
    }
    //classify
    if(!params.skip_illumina_kraken2){
    //if(!params.skip_illumina_kraken2){   
        KRAKEN2_KRAKEN2_ILLUMINA (illumina_reads, PREPARE_REFERENCES.out.ch_kraken2_db, false, true)
        KRAKENTOOLS_COMBINEKREPORTS_ILLUMINA ( KRAKEN2_KRAKEN2_ILLUMINA.out.report.map{ [[id:"kraken2.report"], it[1]] }.groupTuple())
        ch_software_versions = ch_software_versions.mix(KRAKEN2_KRAKEN2_ILLUMINA.out.versions)
    }
    //tools for special organism
    
    SPECIAL_TOOLS_BASED_ON_ILLUMINA(illumina_reads)
    ch_software_versions = ch_software_versions.mix(SPECIAL_TOOLS_BASED_ON_ILLUMINA.out.versions)
    
    

    if(!params.skip_nanopore_reads_qc){
        QC_NANOPORE(nanopore_reads, PREPARE_REFERENCES.out.ch_hostile_ref_minimap2)
        QC_NANOPORE.out.qc_reads
            .filter {meta, reads -> reads.size() > 0 && reads.countFastq() > 0}
            .set { nanopore_reads }
        ch_software_versions = ch_software_versions.mix(QC_NANOPORE.out.versions) 
    }

    //classify
    if(!params.skip_nanopore_kraken2){
        KRAKEN2_KRAKEN2_NANOPORE(
            nanopore_reads, 
            PREPARE_REFERENCES.out.ch_kraken2_db, 
            false, 
            true
        )
        KRAKENTOOLS_COMBINEKREPORTS_NANOPORE(
            KRAKEN2_KRAKEN2_NANOPORE.out.report.map{ [[id:"kraken2_nanopore"], it[1]] }.groupTuple()
        )
        ch_software_versions = ch_software_versions.mix(KRAKEN2_KRAKEN2_NANOPORE.out.versions)
    }

    //tools for special organism
    
    SPECIAL_TOOLS_BASED_ON_NANOPORE(nanopore_reads)
    ch_software_versions = ch_software_versions.mix(SPECIAL_TOOLS_BASED_ON_NANOPORE.out.versions)
    
    
    // assembly
    if(!params.skip_nanopore_reads_assembly){
        
        nanopore_reads.join(QC_NANOPORE.out.qc_stats)//.view()
        nanopore_reads.join(QC_NANOPORE.out.qc_stats).map{
           meta, reads, stats -> [meta, reads, stats.splitCsv(header: true, sep:'\t', strip:true)]
        }.map{
            meta, reads, row -> [meta, reads, row.sum_len[0]]
        }//.view()

        nanopore_reads.join(QC_NANOPORE.out.qc_stats).map{
           meta, reads, stats -> [meta, reads, stats.splitCsv(header: true, sep:'\t', strip:true)]
        }.map{
            meta, reads, row -> [meta, reads, row.sum_len[0]]
        }.branch{
             pass_for_assembly: it[2].toBigInteger() >= params.min_tbp_for_assembly_nanopore
             fail_for_assembly:  it[2].toBigInteger() < params.min_tbp_for_assembly_nanopore
         }.set{  
            ch_input
        }
        
        CSVTK_CONCAT_STATS_NOT_ASSEMBLED(
            ch_input.fail_for_assembly.map{
                meta, reads, row -> [meta, reads]
            }.join(QC_NANOPORE.out.qc_stats).map{
                cfg, reads, stats -> stats
            }.collect().map { files -> tuple([id: "samples_not_assembled"], files)}, 
            in_format, 
            out_format 
        )
        
        nanopore_reads = ch_input.pass_for_assembly.map{
            meta, reads, row -> [meta, reads]
        }
        //nanopore_reads.view()


        //flye with 4x iteration and 1x medaka
        ASSEMBLE_NANOPORE ( nanopore_reads)
        ASSEMBLE_NANOPORE.out.contigs
                //sometime, the file has no content but the size is not zero
                //.filter { meta, contigs -> contigs.size() > 0 }
                .filter { meta, contigs -> contigs.countFasta() > 0 }
                .set { contigs }

        //contigs = ASSEMBLE_NANOPORE.out.contigs
        contig_file_ext = ".fa.gz"
        ch_software_versions = ch_software_versions.mix(ASSEMBLE_NANOPORE.out.versions)
        stats = ASSEMBLE_NANOPORE.out.stats
        
        contigs.view()

        //if(!params.skip_illumina_reads_polish  && params.skip_polypolish && illumina_reads.ifEmpty(null)){
        if(!params.skip_illumina_reads_polish  && !params.skip_polypolish){
               
            
            RUN_POLYPOLISH(illumina_reads, contigs)
            contigs = contigs.merge( RUN_POLYPOLISH.out.contigs.ifEmpty([])).map{
                it -> 
                if(it.size() == 4){
                    [it[2], it[3]]
                }
                else{
                    [it[0], it[1]]
                }
            }.view()
            //print(RUN_POLYPOLISH.out.contigs.getClass())
            //print(RUN_POLYPOLISH.out.contigs.ifEmpty(contigs).getClass())
            
            stats = stats.merge( RUN_POLYPOLISH.out.stats.ifEmpty([])).map{
                it -> 
                if(it.size() == 4){
                    [it[2], it[3]]
                }
                else{
                    [it[0], it[1]]
                }
            }.view()

            ch_software_versions = ch_software_versions.mix(RUN_POLYPOLISH.out.versions)
            contig_file_ext = ".fa.gz"
             
        }
       
        
        if(!params.skip_illumina_reads_polish  && !params.skip_polca){
            
            RUN_POLCA(illumina_reads, contigs)
            //RUN_POLCA.out.contigs.view()
            contigs = contigs.merge( RUN_POLCA.out.contigs.ifEmpty([])).map{
                    it -> 
                    if(it.size() == 4){
                        [it[2], it[3]]
                    }
                    else{
                        [it[0], it[1]]
                    }
                }.view()
            
             stats = stats.merge( RUN_POLCA.out.stats.ifEmpty([])).map{
                    it -> 
                    if(it.size() == 4){
                        [it[2], it[3]]
                    }
                    else{
                        [it[0], it[1]]
                    }
                }.view()

    
            contig_file_ext = ".fa.gz"
            ch_software_versions = ch_software_versions.mix(RUN_POLCA.out.versions)
        }
       
        CSVTK_CONCAT_STATS_ASM(stats.map { cfg, stats -> stats }.collect()
            .map { 
                files -> tuple([id:"assembly_stats"], files)
            }, 
            in_format, 
            out_format 
        ) 
        // analysis

        if(! params.skip_depth_and_coverage){
           
            //contigs.view()
            DEPTH_ILLUMINA(illumina_reads, contigs)
            CSVTK_CONCAT_DEPTH_ILLUMINA(
                DEPTH_ILLUMINA.out.sample_coverage.map { 
                    cfg, stats -> stats }.collect().map { files -> tuple([id:"assembly.depth_illumina"], files)}, 
                in_format, 
                out_format 
            ) 
            DEPTH_NANOPORE(nanopore_reads, contigs)
            CSVTK_CONCAT_DEPTH_NANOPORE(
                DEPTH_NANOPORE.out.sample_coverage.map { cfg, stats -> stats }.collect().map 
                    { files -> tuple([id:"assembly.depth_nanopore"], files)}, 
                in_format, 
                out_format 
            ) 
        }

        if(! params.skip_checkm2){
            ch_input_checkm2 = contigs.map { 
                cfg, contigs -> contigs 
            }.collect().map{
                files -> tuple([id:"checkm2"], files)
            }//.view()
            CHECKM2_PREDICT(ch_input_checkm2, contig_file_ext, PREPARE_REFERENCES.out.ch_checkm2_db) 
        
        }

        
        if(!params.skip_gambit){

            GAMBIT_QUERY(contigs, PREPARE_REFERENCES.out.ch_gambit_db)
            
            ch_input_gambit_query = contigs.map { 
                cfg, contigs -> contigs 
            }.collect().map{
                files -> tuple([id:"gambit_query"], files)
            }

            GAMBIT_QUERY_COLLECT(ch_input_gambit_query, PREPARE_REFERENCES.out.ch_gambit_db)
            
            ch_input_gambit_tree = contigs.map { cfg, contigs -> contigs }.collect()
                 .filter{contigs -> contigs.size() >= 3}
                 .map{files -> tuple([id:"gambit_tree"], files)}//.view()

            GAMBIT_TREE(ch_input_gambit_tree)  
            
        }
        
        ANNOTATION(contigs, PREPARE_REFERENCES.out.ch_bakta_db, PREPARE_REFERENCES.out.ch_amrfinderplus_db)

        SPECIAL_TOOLS_BASED_ON_CONTIGS(contigs)
        ch_software_versions = ch_software_versions.mix(SPECIAL_TOOLS_BASED_ON_CONTIGS.out.versions)
    }
    
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_software_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
