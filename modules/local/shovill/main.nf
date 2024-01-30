process SHOVILL {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shovill:1.1.0--0' :
        'biocontainers/shovill:1.1.0--0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*contigs.fa.gz")                         , emit: contigs
    tuple val(meta), path("shovill.corrections")                , optional: true, emit: corrections
    tuple val(meta), path("shovill.log")                        , emit: log
    tuple val(meta), path("{skesa,spades,megahit,velvet}.fasta"), emit: raw_contigs
    tuple val(meta), path("*contigs.{fastg,gfa,LastGraph}.gz")      , optional:true, emit: gfa
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def memory = task.memory.toGiga()
    """
    shovill \\
        --R1 ${reads[0]} \\
        --R2 ${reads[1]} \\
        $args \\
        --cpus $task.cpus \\
        --ram $memory \\
        --outdir ./ \\
        --force

    mv contigs.fa ${prefix}.contigs.fa
    gzip -n ${prefix}.contigs.fa
    mv contigs.gfa ${prefix}.contigs.gfa
    gzip -n ${prefix}.contigs.gfa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shovill: \$(echo \$(shovill --version 2>&1) | sed 's/^.*shovill //')
    END_VERSIONS
    """
}
