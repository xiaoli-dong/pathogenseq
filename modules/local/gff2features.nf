
process GFF2FEATURES{

    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path('*.tsv'), emit: feature_count

    when:
    task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"

        """
        extract_info_from_gff.pl  -g ${gff} -n ${meta.id} > ${prefix}.feature_count.tsv
        """
   
}
