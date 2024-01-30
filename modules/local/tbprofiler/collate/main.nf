process TBPROFILER_COLLATE {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'
    
    conda "${moduleDir}/environment.yml"
    conda "bioconda::tb-profiler=5.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tb-profiler%3A5.0.1--pyhdfd78af_1' :
        'biocontainers/tb-profiler%3A5.0.1--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(json)

    output:
    tuple val(meta), path("${prefix}.txt") , emit: tsv, optional: true
    tuple val(meta), path("${prefix}.json"), emit: json
    tuple val(meta), path("*.variants.txt") , emit: variants, optional: true
    tuple val(meta), path("*.itol.txt"), emit: itol
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
