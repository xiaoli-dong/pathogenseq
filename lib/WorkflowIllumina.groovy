//
// This file holds several functions specific to the workflow/Illumina.nf in the nf-core/Illumina pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowIllumina {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log, valid_params) {
        WorkflowCommons.genomeExistsError(params, log)
      
        // Generic parameter validation
        if (!valid_params['illumina_reads_assembler'].contains(params.illumina_reads_assembler)) {
            log.error "Invalid option: '${params.illumina_reads_assembler}'. Valid options for '--illumina_reads_assembler': ${valid_params['illumina_read_assembler'].join(', ')}."
            System.exit(1)
        }
        
    }
}
