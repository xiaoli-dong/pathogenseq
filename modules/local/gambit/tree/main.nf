process GAMBIT_TREE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::gambit=1.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gambit%3A1.0.0--py39hbf8eff0_0':
        'biocontainers/gambit%3A1.0.0--py39hbf8eff0_0' }"

    input:
    tuple val(meta), path(contigs)


    output:
    tuple val(meta), path("*.newick") , emit: newick //tree
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''

    def prefix = task.ext.prefix ?: "${meta.id}"


    """
    gambit tree \\
        ${args} \\
        -c ${task.cpus} \\
        ${contigs} > ${prefix}.newick


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gambit: \$(gambit --version 2>&1 | sed -e "s/gambit, version //g")
    END_VERSIONS
    """
}
