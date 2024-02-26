process SKESA {
    tag "$meta.id"
    label 'process_medium'

   
    conda "bioconda::skesa=2.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/skesa:2.5.1--hdcf5f25_0':
        'biocontainers/skesa:2.5.1--hdcf5f25_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.contigs.fa.gz'), emit: contigs
    path ("versions.yml"), emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_reads = meta.single_end ? "-s $reads" : "${reads[0]},${reads[1]}"
    maxmem = task.memory.toGiga()

    """
    skesa $args --reads $input_reads --cores $task.cpus --memory $maxmem --contigs_out ${prefix}.contigs.fa >& ${prefix}_skesa.log
    gzip -m ${prefix}.contigs.fa
    
    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        skesa: \$(skesa -v 2>&1 | sed -e 's/^skesa -v//;;s/^SKESA //;' | sed '/^[[:space:]]*\$/d')
    END_VERSIONS
    
    
    """
    
}
