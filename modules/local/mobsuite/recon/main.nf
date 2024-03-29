
process MOBSUITE_RECON {
    tag "$meta.id"
    label 'process_medium'

   conda "bioconda::mob_suite=3.1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mob_suite:3.1.4--pyhdfd78af_0':
        'biocontainers/mob_suite:3.1.4--pyhdfd78af_0' }"


    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*chromosome.fasta")    , emit: chromosome, optional: true
    tuple val(meta), path("*contig_report.txt")   , emit: contig_report, optional: true
    tuple val(meta), path("*plasmid_*.fasta")     , emit: plasmids        , optional: true
    tuple val(meta), path("*mobtyper_results.txt"), emit: mobtyper_results, optional: true
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    mob_recon \\
        --infile $fasta_name \\
        $args \\
        --num_threads $task.cpus \\
        --outdir results \\
        --prefix $prefix
    
    mv results/* .
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mobsuite: \$(echo \$(mob_recon --version 2>&1) | sed 's/^.*mob_recon //; s/ .*\$//')
    END_VERSIONS
    """
}
