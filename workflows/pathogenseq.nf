/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowPathogenseq.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include {RUN_ILLUMINA_QC} from '../subworkflows/local/qc_illumina'
include {RUN_NANOPORE_QC} from '../subworkflows/local/qc_nanopore'
include {RUN_ASSEMBLE_SHORT} from '../subworkflows/local/assembly_short'
include {
    RUN_ASSEMBLE_LONG;
    //run_dnaA_genome;
} from '../subworkflows/local/assembly_long'

include {
    RUN_POLYPOLISH;
    RUN_POLCA;


} from '../subworkflows/local/short_reads_polisher'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include {
    CSVTK_CONCAT as CONCAT_STATS_SHORT_RAW;
    CSVTK_CONCAT as CONCAT_STATS_SHORT_QC;
    CSVTK_CONCAT as CONCAT_STATS_LONG_RAW;
    CSVTK_CONCAT as CONCAT_STATS_LONG_QC;
    CSVTK_CONCAT as CONCAT_STATS_ASM;
    CSVTK_CONCAT as CONCAT_STATS_BAKTA;
    CSVTK_CONCAT as CONCAT_MOBSUITE;
    CSVTK_CONCAT as CONCAT_AMR;
    CSVTK_CONCAT as CONCAT_MLST;


} from '../modules/nf-core/csvtk/concat/main'

include {KRAKEN2_KRAKEN2} from '../modules/nf-core/kraken2/kraken2/main' 
include {BRACKEN_BRACKEN} from '../modules/nf-core/bracken/bracken/main'
include { BAKTA_BAKTA } from '../modules/nf-core/bakta/bakta/main' 
include {GFF2FEATURES as BAKTA_FEATURES} from '../modules/local/gff2features'  
include {AMRFINDERPLUS_UPDATE} from '../modules/nf-core/amrfinderplus/update/main'
include {AMRFINDERPLUS_RUN} from '../modules/nf-core/amrfinderplus/run/main'
include { MLST } from '../modules/nf-core/mlst/main'
include { MOBSUITE_RECON } from '../modules/nf-core/mobsuite/recon/main'
include { ABRICATE_RUN} from '../modules/nf-core//abricate/run/main'
include { ABRICATE_SUMMARY} from '../modules/nf-core/abricate/summary/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow PATHOGENSEQ {

    ch_software_versions = Channel.empty()
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_software_versions = ch_software_versions.mix(INPUT_CHECK.out.versions)

   

    reads = INPUT_CHECK.out.reads
    short_reads = INPUT_CHECK.out.shortreads
    long_reads = INPUT_CHECK.out.longreads
    in_format = "tsv"
    out_format = "tsv"

    if(!params.skip_short_reads_qc){
        RUN_ILLUMINA_QC(short_reads)
        ch_software_versions = ch_software_versions.mix(RUN_ILLUMINA_QC.out.versions)
        short_reads = RUN_ILLUMINA_QC.out.qc_reads
        CONCAT_STATS_SHORT_RAW(RUN_ILLUMINA_QC.out.raw_stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "short_reads_raw_seqstats"], files)}, in_format, out_format )
        CONCAT_STATS_SHORT_QC(RUN_ILLUMINA_QC.out.qc_stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id: "short_reads_qc_seqstats"], files)}, in_format, out_format )
    }

    if(!params.skip_long_reads_qc){
        RUN_NANOPORE_QC(long_reads)
        long_reads = RUN_NANOPORE_QC.out.qc_reads
        ch_software_versions = ch_software_versions.mix(RUN_NANOPORE_QC.out.versions)
        
        CONCAT_STATS_LONG_RAW(RUN_NANOPORE_QC.out.raw_stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"long_reads_raw_seqstats"], files)}, in_format, out_format )
        CONCAT_STATS_LONG_QC(RUN_NANOPORE_QC.out.qc_stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"long_reads_qc_seqstats"], files)}, in_format, out_format ) 
    }
    //classify
    if(!params.skip_kraken2){
        Channel
            .value(file( "${params.kraken2_db}" ))
            .set { ch_kraken2_db_file }
        KRAKEN2_KRAKEN2 (short_reads, ch_kraken2_db_file, true, true)
        ch_software_versions = ch_software_versions.mix(KRAKEN2_KRAKEN2.out.versions)
        BRACKEN_BRACKEN(KRAKEN2_KRAKEN2.out.report, ch_kraken2_db_file)
        ch_software_versions = ch_software_versions.mix(BRACKEN_BRACKEN.out.versions)
    }

    // assembly
    if(!params.skip_short_reads_assembly && params.assembly_type == 'short'){
    
        RUN_ASSEMBLE_SHORT ( short_reads)
        contigs = RUN_ASSEMBLE_SHORT.out.contigs
        ch_software_versions = ch_software_versions.mix(RUN_ASSEMBLE_SHORT.out.versions)
        stats = RUN_ASSEMBLE_SHORT.out.stats
          
    }

    if(!params.skip_long_reads_assembly && params.assembly_type == 'long'){
        
        //flye with 4x iteration and 1x medaka
        RUN_ASSEMBLE_LONG ( long_reads, short_reads)
        contigs = RUN_ASSEMBLE_LONG.out.contigs
        ch_software_versions = ch_software_versions.mix(RUN_ASSEMBLE_LONG.out.versions)
        stats = RUN_ASSEMBLE_LONG.out.stats
       
        if(!params.skip_short_reads_polish  && !params.skip_polypolish){
            // 4x iterations are recommended
            RUN_POLYPOLISH(short_reads, contigs)
            contigs = RUN_POLYPOLISH.out.assembly
            stats = RUN_POLYPOLISH.out.stats
            ch_software_versions = ch_software_versions.mix(RUN_POLYPOLISH.out.versions)
        }
        // cannot use gzip contig file as input, otherwise it will hang and never finish
        if(!params.skip_short_reads_polish  && !params.skip_polca){

            RUN_POLCA(short_reads, contigs)
            contigs = RUN_POLCA.out.assembly
            stats = RUN_POLCA.out.stats
            ch_software_versions = ch_software_versions.mix(RUN_POLCA.out.versions)
        }

        //ignore this for now becasue 1) the dnaaplier report error when the input has no dnaA gene 
        //run_dnaA_genome(contigs)
        //run_dnaA_genome.out.fasta.view()
    }

    CONCAT_STATS_ASM(stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"assembly_stats"], files)}, in_format, out_format ) 
    // analysis
    
     if(!params.skip_bakta && params.annotation_tool== "bakta"){
        Channel
        .value(file( "${params.bakta_db}" ))
        .set { ch_bakta_db_file }

        BAKTA_BAKTA(contigs, ch_bakta_db_file, [], [])
        BAKTA_FEATURES(BAKTA_BAKTA.out.gff)
        ch_software_versions = ch_software_versions.mix(BAKTA_BAKTA.out.versions)
        CONCAT_STATS_BAKTA(BAKTA_FEATURES.out.feature_count.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"bakta"], files)}, in_format, out_format ) 
        gff = BAKTA_BAKTA.out.gff 
        ffn = BAKTA_BAKTA.out.ffn
        faa = BAKTA_BAKTA.out.faa
    }
    
    //ARG(contigs, ffn, faa)
    if(!params.skip_amr ){
        AMRFINDERPLUS_UPDATE()
        AMRFINDERPLUS_RUN(contigs, AMRFINDERPLUS_UPDATE.out.db)
        ch_software_versions = ch_software_versions.mix(AMRFINDERPLUS_RUN.out.versions)
        CONCAT_AMR(AMRFINDERPLUS_RUN.out.report.map { cfg, amr -> amr }.collect().map { files -> tuple([id:"amrfinderplus"], files)}, in_format, out_format ) 
    }
    if(!params.skip_mlst){
        MLST (contigs)
        ch_software_versions = ch_software_versions.mix(MLST.out.versions)
        CONCAT_MLST (MLST.out.tsv.map { cfg, mlst -> mlst }.collect().map { files -> tuple([id:"mlst"], files)}, in_format, out_format ) 
    }
    if(!params.skip_mobsuite){
        MOBSUITE_RECON (contigs )
        ch_software_versions = ch_software_versions.mix(MOBSUITE_RECON.out.versions)
        CONCAT_MOBSUITE (MOBSUITE_RECON.out.mobtyper_results.map { cfg, plasmid -> plasmid }.collect().map { files -> tuple([id:"mobsuite"], files)}, in_format, out_format ) 
    } 

    if(!params.skip_virulome){
        //virulome
        ABRICATE_RUN(contigs)
        ch_software_versions = ch_software_versions.mix(ABRICATE_RUN.out.versions)
        ABRICATE_SUMMARY ( 
        ABRICATE_RUN.out.report.collect{ meta, report -> report }.map{ report -> [[ id: "virulome"], report]}
        )
    }
    

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_software_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    
   //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowPathogenseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowPathogenseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BRACKEN_BRACKEN.out.reports.collect{it[1]}.ifEmpty([]))
    //ch_multiqc_files.view()

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList() 
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
