process POLYPOLISH {
    tag "$meta.id"
    label 'process_medium'

   
    conda "bioconda::polypolish=0.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/polypolish:0.5.0--hdbdd923_4':
        'biocontainers/polypolish:0.5.0--hdbdd923_4' }"

    input:
    tuple val(meta), path(assembly)
    tuple val(meta), path(sam1), path(sam2)

    output:
    tuple val(meta), path("*.fasta"), emit: assembly
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
   
    """
    polypolish_insert_filter.py \\
        --in1 $sam1 \\
        --in2 $sam2 \\
        --out1 ${prefix}_filtered_1.sam \\
        --out2 ${prefix}_filtered_2.sam

    polypolish \\
        $assembly \\
        ${prefix}_filtered_1.sam \\
        ${prefix}_filtered_2.sam \\
        >${prefix}_polypolish_polished.genome.fasta
    
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polypolish: \$(echo \$(polypolish -V 2>&1) | sed 's/^.*Polypolish //; s/Using.*\$//')
    END_VERSIONS
    """
}
