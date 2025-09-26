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
include { ASSEMBLY_AUTOCYCLER            } from '../subworkflows/local/assembly_autocycler'
include { RUN_POLYPOLISH; RUN_PYPOLCA;    } from '../subworkflows/local/polisher_illumina'
include { DEPTH_ILLUMINA    }       from '../subworkflows/local/depth_illumina'
include { DEPTH_NANOPORE   }       from '../subworkflows/local/depth_nanopore'
include {

    SPECIAL_TOOLS_BASED_ON_CONTIGS;
    SPECIAL_TOOLS_BASED_ON_ILLUMINA;
    SPECIAL_TOOLS_BASED_ON_NANOPORE;
} from '../subworkflows/local/special_tools'

include{GENOMESIZE_NANOPORE} from '../subworkflows/local/genomesize_nanopore'

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
    CSVTK_CONCAT as CSVTK_CONCAT_GTDBTK_ANI_SUMMARY;
    CSVTK_CONCAT as CSVTK_CONCAT_GTDBTK_ANI_CLOSEST;
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
include{
    GTDBTK_ANIREP
} from '../modules/local/gtdbtk/anirep/main'

include {
    GAMBIT_QUERY as GAMBIT_QUERY_COLLECT;
    GAMBIT_QUERY as GAMBIT_QUERY;
 }    from '../modules/local/gambit/query/main'
include { GAMBIT_TREE       }       from '../modules/local/gambit/tree/main'
include {DORADO_VERSION_CHECK} from '../modules/local/dorado/main.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow NANOPORE {

    ch_software_versions = Channel.empty()
    DORADO_VERSION_CHECK()
    //DORADO_VERSION_CHECK.out.versions.view()


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
