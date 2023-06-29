#!/usr/bin/env nextflow

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()
    
include { ABRICATE} from '../../modules/local/abricate' addParams( options: modules['abricate_resistome'])
include {ABRICATE_SUMMARIZE } from '../../modules/local/abricate' addParams( options: modules['abricate_summarize'])


include { SRAX } from '../../modules/local/srax' addParams( options: modules['srax'])
include { RGI} from '../../modules/local/rgi' addParams( options: modules['rgi'])
include { RGI_HEATMAP } from '../../modules/local/rgi' addParams( options: modules['rgi_heatmap'])
include { HAMRONIZE_RGI; HAMRONIZE_SRAX;HAMRONIZE_ABRICATE; HAMRONIZE_SUMMARIZE } from '../../modules/local/hamronization' addParams( options: modules['hamronize'])

//TODO: RGI_HEATMAP is not working, 
workflow ARG {
    take:
        fasta         // channel: [ val(meta), [ fasta ] ]
        card_db

    main:
        //return sampleid_abricate.tsv
        //ABRICATE(fasta)
        //ABRICATE.out.report.collect{ it[1] } | ABRICATE_SUMMARIZE

        SRAX(fasta.collect { it[1] })
        SRAX.out.report | HAMRONIZE_SRAX

        //RGI(fasta, card_db)
        //RGI.out.json.collect{ it[1] } | RGI_HEATMAP

        //ABRICATE.out.report | HAMRONIZE_ABRICATE
        //RGI.out.txt | HAMRONIZE_RGI


       // HAMRONIZE_ABRICATE.out.hamronized.mix(HAMRONIZE_SRAX.out.hamronized).collect{ it[1] } | HAMRONIZE_SUMMARIZE
        HAMRONIZE_SRAX.out.hamronized.collect{ it[1] } | HAMRONIZE_SUMMARIZE
    emit:
        //summary  = ABRICATE_SUMMARIZE.out.summary  // channel: [ summary ]
        out_srax = SRAX.out.report
        //heatmap = RGI_HEATMAP.out.heatmap
        hamronize_summary_tsv =  HAMRONIZE_SUMMARIZE.out.summary_tsv
        hamronize_summary_html =  HAMRONIZE_SUMMARIZE.out.summary_html
}    
