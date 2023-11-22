process SEQKIT_STATS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::seqkit=2.6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit%3A2.6.0--h9ee0642_0':
        'biocontainers/seqkit:2.6.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.tsv"), emit: stats
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gzipped = reads[0].toString().endsWith('.gz')
    def cmd = gzipped ? 'zcat' : 'cat'

    """
    ${cmd} ${reads} | \\
        seqkit stats \\
        -i ${meta.id} \\
        $args \\
        -o ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """
}
