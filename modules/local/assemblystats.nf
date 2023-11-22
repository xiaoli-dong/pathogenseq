process ASSEMBYSTATS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::assembly-stats=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/assembly-stats:1.0.1--h9f5acd7_5':
        'biocontainers/assembly-stats:1.0.1--h9f5acd7_5' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.txt"), emit: stats
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gzipped = fasta.toString().endsWith('.gz')
    def cmd = gzipped ? "<(zcat ${fasta})" : "${fasta}"

    """
    assembly-stats \\
        $args \\
        $cmd \\
        > ${prefix}.txt
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        assembly-stats: \$(echo \$(assembly-stats -v 2>&1) | sed 's/Version: //' ))
    END_VERSIONS
    """

    
}
