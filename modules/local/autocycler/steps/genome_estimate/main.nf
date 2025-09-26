process AUTOCYCLER_GENOMESIZE {
    tag "$meta.id"
    label 'process_medium'

     conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/raven-assembler%3A1.8.3--h5ca1c30_3' :
        'biocontainers/raven-assembler:1.8.3--h5ca1c30_3' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}_genome_size.txt"), emit: genome_size_file
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # run tool
    raven -t $task.cpus $reads > ${prefix}_assembly.fasta

    # Calculate total length of all sequences
    awk '/^>/ {next} {total += length(\$0)} END {print total}' ${prefix}_assembly.fasta > ${prefix}_genome_size.txt


    # get tool version
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        raven: \$( raven --version )
    END_VERSIONS
    """
}

