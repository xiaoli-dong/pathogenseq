process RESTARTGENOME {
    tag "$meta.id"
    label 'process_single'

   
    conda "conda-forge::perl=5.26.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2':
        'biocontainers/perl:5.26.2' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta), path(assembly_info)

    output:
    tuple val(meta), path("*.restart.fasta"), emit: fasta
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    reset_start_position_for_circular_genome.pl \\
        -f $fasta \\
        -i $assembly_info \\
        > $prefix\.restart.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl -v 2>&1) | sed -n \'2 p\' | sed 's/^.*?(v//g; s/^).*//g;' ))
    END_VERSIONS
    """
}
