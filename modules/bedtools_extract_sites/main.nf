process BEDTOOLS_EXTRACT_SITES {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(bedgraph)
    path bed

    output:
    tuple val(meta), path("*.bedGraph"), emit: bedgraph

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedtools intersect \\
        -a ${bedgraph} \\
        -b ${bed} \\
        -wa \\
        -wb \\
        > ${prefix}.bedGraph
    """
}
