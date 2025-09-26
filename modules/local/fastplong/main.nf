process FASTPLONG {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastplong%3A0.3.0--h224cc79_0' :
        'biocontainers/fastplong:0.3.0--h224cc79_0' }"

    input:
    tuple val(meta), path(reads)
    path  adapter_fasta
    val   save_trimmed_fail

    output:
    tuple val(meta), path("${prefix}.fastq.gz"),    emit: reads
    tuple val(meta), path("*.json")           ,     emit: json
    tuple val(meta), path("*.html")           ,     emit: html
    tuple val(meta), path("*.log")            ,     emit: log
    path "versions.yml"                       ,     emit: versions
    tuple val(meta), path("*.fail.fastq.gz")  ,     optional:true, emit: reads_fail

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    adapter_list = adapter_fasta ? "--adapter_fasta ${adapter_fasta}" : ""
    fail_fastq = save_trimmed_fail ? "--failed_out ${prefix}.fail.fastq.gz" : ''
    // Added soft-links to original fastqs for consistent naming in MultiQC
    // Use single ended for interleaved. Add --interleaved_in in config.
   
    """
    
    fastplong \\
        --in ${reads} \\
        --out  ${prefix}.fastq.gz \\
        --thread $task.cpus \\
        --json ${prefix}.json \\
        --html ${prefix}.html \\
        $adapter_list \\
        $fail_fastq \\
        $args \\
        2> >(tee ${prefix}.log >&2)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastplong: \$(fastplong --version 2>&1 | sed -e "s/fastplong //g")
    END_VERSIONS
    """
    
    stub:
   
   
    """
    touch "${prefix}.fastplong.fastq.gz"
    touch "${prefix}.fastplong.json"
    touch "${prefix}.fastplong.html"
    touch "${prefix}.fastplong.log"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastplong: \$(fastplong --version 2>&1 | sed -e "s/fastplong //g")
    END_VERSIONS
    """
}
