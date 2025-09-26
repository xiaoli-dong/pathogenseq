process HOSTILE {
    tag "$meta.id"
    label 'process_medium'
    conda "bioconda::hostile=2.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hostile:2.0.2--pyhdfd78af_0':
        'biocontainers/hostile:2.0.2--pyhdfd78af_0' }"


    

    input:
    tuple val(meta), path(reads)
    val(aligner)
    path(ref) //directory, index has same prefix as the directory name

    output:
    tuple val(meta), path('*clean*fastq.gz'), emit: reads
    tuple val(meta), path('*.log')     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def single_end = reads.flatten().size() == 1 ? true : false
    def input = single_end ? "--fastq1 ${reads[0]}" : "--fastq1 ${reads[0]} --fastq2 ${reads[1]}"
    def simplename = reads[0].getSimpleName()
    
    //print(reads.flatten().size())
    //print(simplename)
    
    if (aligner =~ /bowtie2/) {
        """
        hostile clean \\
            --aligner ${aligner} \\
            --output ./ \\
            --index ${ref}/${ref} \\
            --threads $task.cpus \\
            ${input} \\
            $args \\
            >& ${prefix}.hostile_${aligner}.log

        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hostile: \$(hostile --version)
        END_VERSIONS
        """
    } else if (aligner =~ /minimap2/) {
        
        """
        hostile clean \\
            --aligner ${aligner} \\
            --output ./ \\
            --index ${ref} \\
            --threads $task.cpus \\
            ${input} \\
            $args \\
            >& ${prefix}.hostile_${aligner}.log
        

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hostile: \$(hostile --version)
        END_VERSIONS
        """
    }
}
