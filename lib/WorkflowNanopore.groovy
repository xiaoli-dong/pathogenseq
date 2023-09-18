//
// This file holds several functions specific to the workflow/Nanopore.nf in the nf-core/Nanopore pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowNanopore {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log, valid_params) {
        WorkflowCommons.genomeExistsError(params, log)
        
        // Generic parameter validation
        if (!valid_params['nanopore_reads_assembler'].contains(params.nanopore_reads_assembler)) {
            log.error "Invalid option: '${params.nanopore_read_assember}'. Valid options for '--nanopore_reads_assembler': ${valid_params['nanopore_reads_assembler'].join(', ')}."
            System.exit(1)
        }
    }

    
}
