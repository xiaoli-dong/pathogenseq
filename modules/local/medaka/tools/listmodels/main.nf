process MEDAKA_TOOLS_LISTMODELS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:2.1.0--py38ha0c3a46_0' :
        'biocontainers/medaka:2.1.0--py38ha0c3a46_0' }"

    
    output:
    path("medaka_tools_list_models.txt"), emit: models
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def m_option = mode != "NA" ? "-m ${mode}" : ''
    """
    medaka tools list_models > medaka_tools_list_models.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """
}
