process GBS_SBG {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::blast=2.15.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast%3A2.15.0--pl5321h6f7f691_1':
        'biocontainers/blast%3A2.15.0--pl5321h6f7f691_1' }"


    input:
    tuple val(meta), path(fasta) //contigs
    path(ref)

    output:
    //# Name  Serotype        Uncertainty
    //S18     NT      MaxCov:0;MaxID:0
    //S17     GBS-SBG:Ia
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gzipped = fasta.toString().endsWith('.gz')
    def cmd_input = gzipped ? "zcat ${fasta}" : "cat ${fasta}"
    def cmd_refdb = ref ? "-ref ${ref}" : "" 
    
    """
    ${cmd_input} | GBS-SBG.pl \\
        ${args} \\
        -name ${meta.id} \\
        ${cmd_refdb} \\
        > ${prefix}.tsv
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
