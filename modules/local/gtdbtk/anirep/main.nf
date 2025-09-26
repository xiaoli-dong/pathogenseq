process GTDBTK_ANIREP {
    tag "${meta.id}"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gtdbtk:2.4.1--pyhdfd78af_1':
        'biocontainers/gtdbtk:2.4.1--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(contigs, stageAs: "input_dir/*")
    val fasta_ext
    path(db)


    output:
    //tuple val(meta), path("${prefix}")                               , emit: ani_rep_outdir
    tuple val(meta), path("${prefix}/*.ani_closest.tsv")        , emit: closest
    tuple val(meta), path("${prefix}/*.ani_summary.tsv")      , emit: summary
    tuple val(meta), path("${prefix}/*.json"), emit: json
    tuple val(meta), path("${prefix}/${prefix}.log")                 , emit: log
    tuple val(meta), path("${prefix}/${prefix}.warnings.log")        , emit: warnings
    path ("versions.yml")                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    prefix              = task.ext.prefix ?: "${meta.id}"

    """
    export GTDBTK_DATA_PATH="\$(find -L ${db} -name 'metadata' -type d -exec dirname {} \\;)"
    gtdbtk ani_rep \\
        ${args} \\
        --genome_dir input_dir \\
        --out_dir ${prefix} \\
        --prefix ${meta.id} \\
        -x $fasta_ext \\
        --cpus ${task.cpus}

    mv ${prefix}/gtdbtk.log "${prefix}/${prefix}.log"
    mv ${prefix}/gtdbtk.warnings.log "${prefix}/${prefix}.warnings.log"
    mv ${prefix}/gtdbtk.json "${prefix}/${prefix}.json"


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo \$(gtdbtk --version 2>/dev/null) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}

    touch ${prefix}/${meta.id}.ani_summary.tsv
    touch ${prefix}/${meta.id}.ani_closest.tsv
    touch ${prefix}/${prefix}.log
    touch ${prefix}/${prefix}.json
    touch ${prefix}/${prefix}.warnings.log


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo \$(gtdbtk --version 2>/dev/null) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """
}
