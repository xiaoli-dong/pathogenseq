process TBPROFILER_CUSTOMIZE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::tb-profiler=4.4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tb-profiler:4.4.2--pyh7cba7a3_0 ' :
        'biocontainers/tb-profiler:4.4.2--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("bam/*.bam")     , emit: bam
    tuple val(meta), path("results/*.csv") , emit: csv, optional: true
    tuple val(meta), path("results/*.json"), emit: json
    tuple val(meta), path("results/*.txt") , emit: txt, optional: true
    tuple val(meta), path("vcf/*.vcf.gz")  , emit: vcf
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def input_reads = meta.single_end ? "--read1 $reads" : "--read1 ${reads[0]} --read2 ${reads[1]}"
    def memory = task.memory.toGiga()
    """
    
    tb-profiler \\
        profile \\
        --external_db $db \\
        $args \\
        --prefix ${prefix} \\
        --threads $task.cpus \\
        --ram $memory \\
        --temp ./ \\
        $input_reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tbprofiler:  \$( echo \$(tb-profiler version 2>&1) | sed 's/TBProfiler version //')
    END_VERSIONS
    """
}
