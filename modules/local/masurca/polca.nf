process MASURCA_POLCA {
    tag "$meta.id"
    label 'process_medium'

    
    conda "bioconda::masurca=4.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/masurca:4.1.0--pl5321hb5bd705_1':
        'biocontainers/masurca:4.1.0--pl5321hb5bd705_1' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("*.fa.gz"), emit: contigs
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gzipped = contigs.toString().endsWith('.gz')
    def outfile = gzipped ? file(contigs.baseName).baseName : contigs.baseName
    def command = gzipped ? 'zcat' : 'cat'
    
    // cannot use gzip contig file as  for polca, otherwise it will hang and never finish

    """
    $command $contigs > ${outfile}.fixed.fa

    polca.sh \\
    $args \\
    -a ${outfile}.fixed.fa \\
    -r \'$reads\' \\
    -t ${task.cpus}
    
    mv ${outfile}.fixed.fa.PolcaCorrected.fa ${prefix}.fa
    gzip -n ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        masurca: \$(echo \$(masurca --version 2>&1) | sed 's/^.*version //;')
    END_VERSIONS
    """
}
