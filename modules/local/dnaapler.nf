process DNAAPLER {
    tag "$meta.id"
    label 'process_single'
    errorStrategy 'ignore'

    conda "bioconda::dnaapler=0.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dnaapler:0.1.0--pyhdfd78af_0':
        'biocontainers/dnaapler:0.1.0--pyhdfd78af_0' }"

    input:
    
    tuple val(meta), val(seqid), path(fasta)

    output:
    
    tuple val(meta), val(seqid), path("*_reoriented.fasta"), emit: fasta
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
   
    
    """
    dnaapler \\
        chromosome \\
        -i $fasta \\
        #-o . \\
        #-f \\
        -p $prefix\_${seqid} 
    
   
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dnaapler: \$(echo \$(dnaapler --version 2>&1) | sed 's/^.*version //;' ))
    END_VERSIONS
    """
}
