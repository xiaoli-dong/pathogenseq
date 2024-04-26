/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def valid_params = [
    illumina_reads_assembler : ['unicycler', 'spades', 'skesa','megahit', 'shovill']
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
 
//print(params)
// Validate input parameters
WorkflowIllumina.initialise(params, log, valid_params)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
//def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]

//for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include {
    INPUT_CHECK 
} from '../subworkflows/local/input_check'
include {
    QC_ILLUMINA
} from '../subworkflows/local/qc_illumina'
include {
    ASSEMBLE_ILLUMINA
} from '../subworkflows/local/assembly_illumina'
include {
    ANNOTATION
} from '../subworkflows/local/annotation'
include {
    PREPARE_REFERENCES
} from '../subworkflows/local/prepare_references'
include {
    DEPTH_ILLUMINA
}       from '../subworkflows/local/depth_illumina'

include { 
    SPECIAL_TOOLS_BASED_ON_ILLUMINA;
    SPECIAL_TOOLS_BASED_ON_CONTIGS;
} from '../subworkflows/local/special_tools'

/*
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include {
    CUSTOM_DUMPSOFTWAREVERSIONS
} from '../modules/nf-core/custom/dumpsoftwareversions/main'
include {
    KRAKEN2_KRAKEN2 as KRAKEN2_KRAKEN2_ILLUMINA;
} from '../modules/nf-core/kraken2/kraken2/main' 
include {
    KRAKENTOOLS_COMBINEKREPORTS as KRAKENTOOLS_COMBINEKREPORTS_ILLUMINA;
} from '../modules/nf-core/krakentools/combinekreports/main'
include {
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_ASM; 
    CSVTK_CONCAT as CSVTK_CONCAT_DEPTH_ILLUMINA;    
    CSVTK_CONCAT as CSVTK_CONCAT_STATS_NOT_ASSEMBLED;
   
} from '../modules/nf-core/csvtk/concat/main'
/* include {
    PNEUMOCAT
} from '../modules/nf-core/pneumocat/main.nf'

 */
//
// MODULE: local modules
//
include { 
    CHECKM2_PREDICT
} from '../modules/local/checkm2/predict.nf'
include { 
    GAMBIT_QUERY as GAMBIT_QUERY_COLLECT;
    GAMBIT_QUERY as GAMBIT_QUERY;
 } from '../modules/local/gambit/query/main'
include { 
    GAMBIT_TREE
} from '../modules/local/gambit/tree/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow ILLUMINA {

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
    //reads.view()
    illumina_reads = INPUT_CHECK.out.shortreads
    //illumina_reads.view()
    in_format = "tsv"
    out_format = "tsv"
    contig_file_ext = ".fa.gz"
    
    //illumina_reads.view()
    if(!params.skip_illumina_reads_qc){

        QC_ILLUMINA(
            illumina_reads,
            [],
            PREPARE_REFERENCES.out.ch_hostile_ref_bowtie2
        )
        ch_software_versions = ch_software_versions.mix(QC_ILLUMINA.out.versions)
        //QC_ILLUMINA.out.qc_reads.view()
        
        //get rid of zero size contig file and avoid the downstream crash
        QC_ILLUMINA.out.qc_reads
            .filter {meta, reads -> reads[0].size() > 0 && reads[0].countFastq() > 0}
            .set { illumina_reads }
    }

    //classify
    if(!params.skip_illumina_kraken2){
        KRAKEN2_KRAKEN2_ILLUMINA(
            illumina_reads, 
            PREPARE_REFERENCES.out.ch_kraken2_db, 
            false, 
            true
        )
        KRAKENTOOLS_COMBINEKREPORTS_ILLUMINA(
            KRAKEN2_KRAKEN2_ILLUMINA.out.report.map{ [[id:"kraken2_illumina"], it[1]] }.groupTuple()
        )
        ch_software_versions = ch_software_versions.mix(KRAKEN2_KRAKEN2_ILLUMINA.out.versions)
    }

    //tools for special organism
    SPECIAL_TOOLS_BASED_ON_ILLUMINA(illumina_reads)
    ch_software_versions = ch_software_versions.mix(SPECIAL_TOOLS_BASED_ON_ILLUMINA.out.versions)

    // assembly
    if(!params.skip_illumina_reads_assembly){
        //illumina_reads.view()
        
        illumina_reads.join(QC_ILLUMINA.out.qc_stats).map{
           meta, reads, stats -> [meta, reads, stats.splitCsv(header: true, sep:'\t', strip:true)]
        }.map{
            meta, reads, row -> [meta, reads, row.sum_len[0]]
        }.branch{
             pass_for_assembly: it[2].toBigInteger() >= params.min_tbp_for_assembly_illumina
             fail_for_assembly:  it[2].toBigInteger() < params.min_tbp_for_assembly_illumina
         }.set{  
            ch_input
        }
       
        CSVTK_CONCAT_STATS_NOT_ASSEMBLED(
            ch_input.fail_for_assembly.map{
                meta, reads, row -> [meta, reads]
            }.join(QC_ILLUMINA.out.qc_stats).map{
                cfg, reads, stats -> stats
            }.collect().map { files -> tuple([id: "samples_not_assembled"], files)}, 
            in_format, 
            out_format 
        )
        illumina_reads = ch_input.pass_for_assembly.map{
            meta, reads, row -> [meta, reads]
        }
        //illumina_reads.view()

        ASSEMBLE_ILLUMINA ( illumina_reads)
        // zero size contig can cause some of the program such as bakta, mobsuite file
        ASSEMBLE_ILLUMINA.out.contigs
                .filter { meta, contigs -> contigs.countFasta() > 0 }
                .set { contigs }
        contig_file_ext = ASSEMBLE_ILLUMINA.out.contig_file_ext
        ch_software_versions = ch_software_versions.mix(ASSEMBLE_ILLUMINA.out.versions)
        stats = ASSEMBLE_ILLUMINA.out.stats
        CSVTK_CONCAT_STATS_ASM(
            stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"assembly_stats"], files)}, 
            in_format, 
            out_format
        ) 
       
        illumina_reads.join(contigs)//.view()
        if(! params.skip_depth_and_coverage){
            illumina_reads.join(contigs).multiMap{
                it ->
                reads: [it[0], it[1]]
                contigs: [it[0], it[2]]
            }.set{
                ch_input_depth
            }
            DEPTH_ILLUMINA(ch_input_depth)
            CSVTK_CONCAT_DEPTH_ILLUMINA(
                DEPTH_ILLUMINA.out.sample_coverage.map { 
                    cfg, stats -> stats 
                }.collect().map { 
                    files -> tuple([id:"assembly.depth_illumina"], files)
                    }, 
                in_format, 
                out_format 
            ) 
        }

        if(! params.skip_checkm2){
            ch_input_checkm2 = contigs.map { cfg, contigs -> contigs }.collect().map{files -> tuple([id:"checkm2"], files)}//.view()
            CHECKM2_PREDICT(
                ch_input_checkm2, 
                contig_file_ext, 
                PREPARE_REFERENCES.out.ch_checkm2_db
            ) 
            ch_software_versions = ch_software_versions.mix(CHECKM2_PREDICT.out.versions)
        }

        if(! params.skip_gambit){
            GAMBIT_QUERY(contigs, PREPARE_REFERENCES.out.ch_gambit_db)
            GAMBIT_QUERY_COLLECT(
                contigs.map { cfg, contigs -> contigs }.collect().map{files -> tuple([id:"gambit_query"], files)},
                PREPARE_REFERENCES.out.ch_gambit_db
            )
            
            GAMBIT_TREE(
                contigs.map { 
                    cfg, contigs -> contigs 
                }.collect().filter{contigs -> contigs.size() >= 3}.map{files -> tuple([id:"gambit_tree"], files)}
            )  
        }

        ANNOTATION(
            contigs, 
            PREPARE_REFERENCES.out.ch_bakta_db, 
            PREPARE_REFERENCES.out.ch_amrfinderplus_db
        )
        ch_software_versions = ch_software_versions.mix(ANNOTATION.out.versions)

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
