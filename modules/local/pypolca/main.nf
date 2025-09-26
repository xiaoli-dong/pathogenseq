process PYPOLCA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypolca:0.3.1--pyhdfd78af_1':
        'biocontainers/pypolca:0.3.1--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("out_dir/*_corrected.fasta.gz"), emit: corrected_contigs
    tuple val(meta), path("out_dir/*.report"), emit: report
    tuple val(meta), path("out_dir/*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def gzipped = contigs.toString().endsWith('.gz')
    def reads_cmd = meta.single_end ? "$reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    def outfile = gzipped ? file(contigs.baseName).baseName : contigs.baseName
    def command = gzipped ? 'zcat' : 'cat'

    // cannot use gzip contig file as  for polca, otherwise it will hang and never finish

    """
    $command $contigs > ${outfile}.fixed.fa
    pypolca run --careful \\
        $args \\
        -a  ${outfile}.fixed.fa ${reads_cmd} \\
        -t ${task.cpus} -p ${prefix} \\
        -o out_dir


    gzip -n out_dir/${prefix}_corrected.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypolca: \$(echo \$(pypolca -V 2>&1) | sed 's/^.*version //;')
    END_VERSIONS
    """
}
