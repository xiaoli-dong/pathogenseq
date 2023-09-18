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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_STATS_ASM; } from '../modules/nf-core/csvtk/concat/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_KRAKEN2_ILLUMINA } from '../modules/nf-core/kraken2/kraken2/main' 
include { KRAKENTOOLS_COMBINEKREPORTS as KRAKENTOOLS_COMBINEKREPORTS_ILLUMINA } from '../modules/nf-core/krakentools/combinekreports/main.nf'

//MODULES: local modules
include {CHECKM2_PREDICT} from '../modules/local/checkm2/predict.nf'


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
        QC_ILLUMINA(illumina_reads)
        ch_software_versions = ch_software_versions.mix(QC_ILLUMINA.out.versions)
        
        //get rid of zero size contig file and avoid the downstream crash
        QC_ILLUMINA.out.qc_reads
            .filter {meta, reads -> reads[0].size() > 0 && reads[0].countFastq() > 0}
            .set { illumina_reads }
    }
    //classify
    if(!params.skip_illumina_kraken2){
        KRAKEN2_KRAKEN2_ILLUMINA (illumina_reads, PREPARE_REFERENCES.out.ch_kraken2_db, false, true)
        KRAKENTOOLS_COMBINEKREPORTS_ILLUMINA ( KRAKEN2_KRAKEN2_ILLUMINA.out.report.map{ [[id:"kraken2.report"], it[1]] }.groupTuple())
        ch_software_versions = ch_software_versions.mix(KRAKEN2_KRAKEN2_ILLUMINA.out.versions)
    }

    if(!params.skip_nanopore_reads_qc){
        QC_NANOPORE(nanopore_reads)
        QC_NANOPORE.out.qc_reads
            .filter {meta, reads -> reads.size() > 0 && reads.countFastq() > 0}
            .set { nanopore_reads }
        ch_software_versions = ch_software_versions.mix(QC_NANOPORE.out.versions) 
    }

    // assembly
    if(!params.skip_nanopore_reads_assembly){
        
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
       
        if(!params.skip_illumina_reads_polish  && !params.skip_polypolish){
            // 4x iterations are recommended
            contigs.view()
            illumina_reads.view()

            illumina_reads.join(contigs).multiMap{
                it ->
                illumina_reads: [it[0], it[1]]
                contigs: [it[0], it[2]]
            }.set{
                ch_input_polypolish
            }

            RUN_POLYPOLISH(ch_input_polypolish.illumina_reads, ch_input_polypolish.contigs)
            //RUN_POLYPOLISH(illumina_reads, contigs)
            RUN_POLYPOLISH.out.contigs
                //.filter { meta, contigs -> contigs.size() > 0 }
                .filter { meta, contigs -> contigs.countFasta() > 0 }
                .set { contigs }

            //contigs = RUN_POLYPOLISH.out.contigs
            stats = RUN_POLYPOLISH.out.stats
            ch_software_versions = ch_software_versions.mix(RUN_POLYPOLISH.out.versions)
            contig_file_ext = ".fa.gz"
        }

        
        if(!params.skip_illumina_reads_polish  && !params.skip_polca){

            illumina_reads.join(contigs).multiMap{
                it ->
                illumina_reads: [it[0], it[1]]
                contigs: [it[0], it[2]]
            }.set{
                ch_input_polca
            }
            RUN_POLCA(ch_input_polca.illumina_reads, ch_input_polca.contigs)
            RUN_POLCA.out.contigs
                //.filter { meta, contigs -> contigs.size() > 0 }
                .filter { meta, contigs -> contigs.countFasta() > 0 }
                .set { contigs }

            //contigs = RUN_POLCA.out.contigs
            stats = RUN_POLCA.out.stats
            contig_file_ext = ".fa.gz"
            ch_software_versions = ch_software_versions.mix(RUN_POLCA.out.versions)
        }
        
       CSVTK_CONCAT_STATS_ASM(stats.map { cfg, stats -> stats }.collect().map { files -> tuple([id:"contigs.seqstats"], files)}, in_format, out_format ) 
        // analysis

        if(! params.skip_checkm2){
            ch_input_checkm2 = contigs.map { cfg, contigs -> contigs }.collect().map{files -> tuple([id:"checkm2"], files)}//.view()
            CHECKM2_PREDICT(ch_input_checkm2, contig_file_ext, PREPARE_REFERENCES.out.ch_checkm2_db) 
        
        }

        ANNOTATION(contigs, PREPARE_REFERENCES.out.ch_bakta_db, PREPARE_REFERENCES.out.ch_amrfinderplus_db)
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
