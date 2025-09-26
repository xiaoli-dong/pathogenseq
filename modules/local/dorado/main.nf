process DORADO_VERSION_CHECK {

    label 'process_high'

    container 'nanoporetech/dorado:latest'  // Always pin specific version

    cpus 1
    output:
        path "versions.txt"
    script:
    """
    dorado --version 2>&1 | head -n1 | sed 's/^/dorado,/' >> versions.txt
    minimap2 --version | head -n 1 | sed 's/^/minimap2,/' >> versions.txt
    """
}

/*  dorado polish requires input alignments to be produced using Dorado aligner. 
    When Dorado aligner outputs alignments to stdout, they are not sorted automatically. 
    Instead, samtools needs to be used to sort and index the BAM file
*/
process DORADO_POLISH {

    label 'process_high'

    container 'nanoporetech/dorado:latest'  // Always pin specific version

    input:
        contigs
    output:
        path "versions.txt"
    script:
    """
    dorado --version 2>&1 | head -n1 | sed 's/^/dorado,/' >> versions.txt
    minimap2 --version | head -n 1 | sed 's/^/minimap2,/' >> versions.txt
    """
}

process DORADO_ALIGNER {

    label 'process_high'

    container 'nanoporetech/dorado:latest'  // Always pin specific version

    input:
        tuple val(meta), path(contigs)
        tuple val(meta), path(read_bam)

    output:
        path "versions.txt"

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    
    dorado aligner ${contigs} ${reads.bam} | samtools sort --threads $task.cpus > ${prefix}.dorado_aln.bam
    samtools index ${prefix}.dorado_aln.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -n '/^[0-9]/p')
    END_VERSIONS
    """
}

