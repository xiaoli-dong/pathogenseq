process CHECKM2_PREDICT {
    tag "$meta.id"
    label 'process_single'

    
    conda "bioconda::checkm2=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkm2:1.0.1--pyh7cba7a3_0':
        'biocontainers/checkm2:1.0.1--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(contigs, stageAs: "input_dir/*")
    val fasta_ext
    path db

    output:
    tuple val(meta), path("*.tsv"), emit: report
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when
   
    //checkm2 predict --input input -o output --database_path /nfs/APL_Genomics/db/prod/CheckM2_database/uniref100.KO.1.dmnd  --force -x .fa.gz --allmodels -t 8
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    checkm2 \\
        predict \\
        $args \\
        --input input_dir \\
        -x $fasta_ext \\
        --output-directory ./output \\
        --database_path $db \\
        -t $task.cpus 
    mv ./output/quality_report.tsv ${prefix}.quality_report.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(echo \$(checkm2 --version 2>&1))
    END_VERSIONS
    """

}
