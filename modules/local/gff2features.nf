process GFF2FEATURES{

    tag "$meta.id"
    label 'process_medium'

   
    conda "conda-forge::perl=5.26.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2':
        'biocontainers/perl:5.26.2' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path('*.tsv'), emit: tsv //feature_count
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    extract_info_from_gff.pl  -g ${gff} -n ${meta.id} > ${prefix}.feature_count.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl -v 2>&1) | sed -n \'2 p\' | sed 's/^.*?(v//g; s/^).*//g;' ))
    END_VERSIONS
    """
   
}
