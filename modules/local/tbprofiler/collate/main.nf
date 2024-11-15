process TBPROFILER_COLLATE {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'
    
    conda "${moduleDir}/environment.yml"
    conda "bioconda::tb-profiler=6.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tb-profiler%3A6.2.1--pyhdfd78af_0' :
        'biocontainers/tb-profiler%3A6.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(json)


    output:
    tuple val(meta), path("${prefix}.txt") , emit: txt
    //tuple val(meta), path("${prefix}.json"), emit: json
    tuple val(meta), path("${prefix}.variants.txt") , emit: variants_txt
    tuple val(meta), path("*.itol.txt"), emit: itol, optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    
    
    """
    tb-profiler \\
        collate \\
        $args \\
        --prefix ${prefix} \\
        --dir ./
        

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tbprofiler:  \$( echo \$(tb-profiler --version 2>&1) | sed 's/TBProfiler version //')
    END_VERSIONS
    """
}
