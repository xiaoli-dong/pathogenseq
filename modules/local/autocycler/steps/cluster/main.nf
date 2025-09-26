// modules/autocycler/cluster.nf
process AUTOCYCLER_CLUSTER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/autocycler%3A0.5.0--h3ab6199_0' :
        'autocycler:0.5.0--h3ab6199_0' }"


    input:
     tuple val(meta), path (autocycler_dir)

    output:
    tuple val(meta), path("${autocycler_dir}/clustering"), emit: cluster

    script:
    """
    autocycler cluster -a ${autocycler_dir} --max_contigs 10000 2>> autocycler_cluster.stderr
    """
}
