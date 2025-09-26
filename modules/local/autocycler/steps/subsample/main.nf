process AUTOCYCLER_SUBSAMPLE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/autocycler%3A0.5.0--h3ab6199_0' :
        'autocycler:0.5.0--h3ab6199_0' }"

    input:
    tuple val(meta), path(reads), path(genome_size_file)

    output:
    tuple val(meta), path("${prefix}_subsampled_reads"), emit: out
    tuple val(meta), path("${prefix}_subsample_skipped.txt"), optional: true, emit: skipped
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    // Extract the genome size value from the file
    //def genome_size = genome_size_file.text.trim()

    //--min_read_depth ${min_read_depth} --count ${count} \\
    //AUTOCYCLER_SUBSAMPLE.out
   // .filter { meta, path -> path.list().size() > 0 }
    //| run_flye

    """
    genome_size=\$(cat ${genome_size_file} | tr -d '[:space:]')
    autocycler subsample \\
        $args \\
        --reads ${reads} \\
        --out_dir ${prefix}_subsampled_reads \\
        --genome_size \$genome_size \\
        2>> ${prefix}_autocycler.stderr
    status=\$?

    if grep -q "too shallow to subset" ${prefix}_autocycler.stderr || [ \$status -ne 0 ]; then
        echo "Subsampling failed due to shallow coverage for sample: ${meta.id}" > ${prefix}_subsample_skipped.txt
        # Clean up failed output to avoid confusion
        rm -rf ${prefix}_subsampled_reads/*
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$( autocycler -V )
    END_VERSIONS
    """
}
