process COMBINE_XML {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"
    
    input:
    tuple val(meta), path(xmls)
   
    output:
    tuple val(meta), path("${prefix}.xml"), emit: xml
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    print(xmls)
    """
    combine_xml.py \\
        --output ${prefix}.xml \\
        --input '${xmls}'
       
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
