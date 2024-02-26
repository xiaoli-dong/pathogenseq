process POLYPOLISH {
    tag "$meta.id"
    label 'process_medium'

   
    conda "bioconda::polypolish=0.6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/polypolish:0.6.0--hdbdd923_0':
        'biocontainers/polypolish:0.6.0--hdbdd923_0' }"

    input:
    tuple val(meta), path(draft_contigs)
    tuple val(meta), path(sam1), path(sam2)

    output:
    tuple val(meta), path("*.fasta.gz"), emit: contigs
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
   
    """
    polypolish filter \\
        ${args} \\
        --in1 ${sam1} \\
        --in2 ${sam2} \\
        --out1 ${prefix}.filtered.R1.sam \\
        --out2 ${prefix}.filtered.R2.sam

    polypolish polish \\
        ${args2} \\
        ${draft_contigs} \\
        ${prefix}.filtered.R1.sam \\
        ${prefix}.filtered.R2.sam \\
        > ${prefix}.fasta

    gzip -n ${prefix}.fasta
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polypolish: \$(echo \$(polypolish -V 2>&1) | sed 's/^.*Polypolish //; s/Using.*\$//')
    END_VERSIONS
    """
}
