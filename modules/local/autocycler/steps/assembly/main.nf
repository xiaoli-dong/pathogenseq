process AUTOCYCLER_ASSEMBLY {

    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/autocycler%3A0.5.0--h3ab6199_0' :
        'autocycler:0.5.0--h3ab6199_0' }"

    input:
    tuple val(meta), path (assemblies_dir)

    output:
    tuple val(meta), path("${prefix}_autocycler_out"), emit: autocycler_out

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"


    """
    max_contigs=10000
    # Step 3: compress the input assemblies into a unitig graph
    autocycler compress --max_contigs \$max_contigs -i ${assemblies_dir} -a ${prefix}_autocycler_out2>> autocycler.stderr

    # Capture the exit code
    compress_exit_code=\$?
    # Check success or failure
    if [ \$compress_exit_code -ne 0 ]; then
        echo "Error: autocycler compress failed with exit code \$compress_exit_code" >&2
        exit \$compress_exit_code
    fi

    # Step 4: cluster the input contigs into putative genomic sequences
    autocycler cluster --max_contigs \$max_contigs -a ${prefix}_autocycler_out2>> autocycler.stderr

    if [ ! -d "${prefix}_autocycler_out/clustering/qc_pass/" ]; then
        echo "${prefix}_autocycler_out/clustering/qc_pass/ directory does not exist."
        exit
    fi
    # Steps 5 and 6: trim and resolve each QC-pass cluster
    for c in ${prefix}_autocycler_out/clustering/qc_pass/cluster_*; do
        autocycler trim -c \$c 2>> autocycler.stderr
        autocycler resolve -c \$c 2>> autocycler.stderr
    done

    # Step 7: combine resolved clusters into a final assembly
    autocycler combine -a ${prefix}_autocycler_out -i ${prefix}_autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa 2>> autocycler.stderr


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$( autocycler -V )
    END_VERSIONS

    """
}
