process HOSTILE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::hostile=0.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hostile%3A0.4.0--pyhdfd78af_0':
        'biocontainers/hostile%3A0.4.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    val(aligner)
    path(ref) //directory, index has same prefix as the directory name

    output:
    tuple val(meta), path('*.dehost*fastq.gz'), emit: reads
    tuple val(meta), path('*.log')     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    //def input = meta.single_end ? "--fastq1 ${reads[0]}" : "--fastq1 ${reads[0]} --fastq2 ${reads[1]}"
    def single_end = reads.flatten().size() == 1 ? true : false
    def input = single_end ? "--fastq1 ${reads[0]}" : "--fastq1 ${reads[0]} --fastq2 ${reads[1]}"
    def simplename = reads[0].getSimpleName()
    
    //print(reads.flatten().size())
    //print(simplename)
    
    if (aligner =~ /bowtie2/) {
        """
        hostile clean \\
            --aligner ${aligner} \\
            --out-dir ./ \\
            --index ${ref}/${ref} \\
            --threads $task.cpus \\
            ${input} \\
            $args \\
            >& ${prefix}.hostile_${aligner}.log

        if [[ ${single_end} = true ]]
        then
            mv ${simplename}*.clean.fastq.gz ${prefix}.dehost.fastq.gz
        else
            mv ${simplename}*.clean_1.fastq.gz ${prefix}.dehost.R1.fastq.gz
            mv ${simplename}*.clean_2.fastq.gz ${prefix}.dehost.R2.fastq.gz
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hostile: \$(hostile --version)
        END_VERSIONS
        """
    } else if (aligner =~ /minimap2/) {
        
        """
        hostile clean \\
            --aligner ${aligner} \\
            --out-dir ./ \\
            --index ${ref} \\
            --threads $task.cpus \\
            ${input} \\
            $args \\
            >& ${prefix}.hostile_${aligner}.log
        
        if [[ ${single_end} = true ]]
        then
            mv ${simplename}*.clean.fastq.gz ${prefix}.dehost.fastq.gz
        else
            mv ${simplename}*.clean_1.fastq.gz ${prefix}.dehost.R1.fastq.gz
            mv ${simplename}*.clean_2.fastq.gz ${prefix}.dehost.R2.fastq.gz
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hostile: \$(hostile --version)
        END_VERSIONS
        """
    }
}
