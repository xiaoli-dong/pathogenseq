process FORMATCSV {
    tag "$meta.id"
    label 'process_low'

    conda 'bioconda::shiptv=0.4.0'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/shiptv:0.4.0--pyh5e36f6f_0'
    } else {
    container 'quay.io/biocontainers/shiptv:0.4.0--pyh5e36f6f_0'
    }

    input:
    tuple val(meta), path(csv)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    formatCSV.py $args -n ${meta.id} -o ${prefix}.tsv $csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
