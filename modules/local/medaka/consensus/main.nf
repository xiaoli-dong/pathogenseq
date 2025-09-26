process MEDAKA_CONSENSUS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:2.1.0--py38ha0c3a46_0' :
        'biocontainers/medaka:2.1.0--py38ha0c3a46_0' }"

    input:
    tuple val(meta), path(reads), path(assembly)
    val mode

    output:
    tuple val(meta), path("medaka/*.fasta.gz"), emit: assembly
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def m_option = mode != "NA" ? "-m ${mode}" : ''
    """
    medaka_consensus \\
        -t $task.cpus \\
        $args \\
        ${m_option} \\
        -i $reads \\
        -d $assembly

    mv medaka/consensus.fasta medaka/${prefix}.fasta

    gzip -n medaka/${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """
}
