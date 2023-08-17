process BRACKEN_BRACKEN {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI.
    // Please update version string below when bumping container versions.
    //for some the dataset bracken 2.7 gives "Error: no reads found. Please check your Kraken report"
    
    conda "bioconda::bracken=2.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bracken%3A2.8--py39h1f90b4d_1':
        'biocontainers/bracken:2.8--py39h1f90b4d_1' }"

    input:
    tuple val(meta), path(kraken_report)
    path database

    output:
    tuple val(meta), path(bracken_report), emit: reports
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    bracken_report = "${prefix}.tsv"
    // WARN: Version information not provided by tool on CLI.
    // Please update version string below when bumping container versions.
    def VERSION = '2.7'
    """
    bracken \\
        ${args} \\
        -d '${database}' \\
        -i '${kraken_report}' \\
        -o '${bracken_report}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bracken: ${VERSION}
    END_VERSIONS
    """
}
