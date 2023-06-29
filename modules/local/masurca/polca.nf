process MASURCA_POLCA {
    tag "$meta.id"
    label 'process_high'

    
    conda "bioconda::masurca=4.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/masurca:4.1.0--pl5321hb5bd705_1':
        'biocontainers/masurca:4.1.0--pl5321hb5bd705_1' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.fa"), emit: assembly
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    polca.sh \\
    -a ${assembly} \\
    -r \'$reads\' \\
    -t ${task.cpus}
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(masurca --version 2>&1) | sed 's/^.*version //;))
    END_VERSIONS
    """
}
