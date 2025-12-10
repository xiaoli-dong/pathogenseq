VERSION = '1.2.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
process PNEUMOCAT {
    tag "$meta.id"
    label 'process_single'
    errorStrategy 'ignore'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pneumocat:1.2.1--0':
        'biocontainers/pneumocat:1.2.1--0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.final_results.xml"), emit: results
    tuple val(meta), path("coverage_summary.txt"), emit: coverage
    path "versions.yml"           , emit: versions
    /*
    If only one capsular type is matched with more than 90% coverage 
    then the report from step 1 contained in this xml file is considered 
    the final result (result type="Serotype") and no further folders 
    will appear within the PneumoCaT output folder. If more than one 
    capsular type are matched with more than 90% coverage then the 
    software moves to step two and a SNP_based_serotyping folder is 
    created containing a second XML file with the final result 
    - see STEP 2- VARIANT-BASED APPROACH.
    Note that the output XML file from step 1 only reports two capsular types, 
    when actually more could be matched and all will pass to step 2 for 
    further distinction. Further information on mapped serotypes in 
    stage 1 can be found in "Coverage_summary.txt". If the top hit 
    coverage is < 90% then no serotypes are reported and 'Failed' 
    appears instead.
    */
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    PneumoCaT.py \\
        --input_directory ./ \\
        $args \\
        --threads $task.cpus \\
        --output_dir ./

    if [ -d "SNP_based_serotyping" ]
    then
        cp SNP_based_serotyping/${prefix}.results.xml ${prefix}.final_results.xml
    else
        cp ${prefix}.results.xml ${prefix}.final_results.xml
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pneumocat: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.results.xml
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pneumocat: $VERSION
    END_VERSIONS
    """
}
