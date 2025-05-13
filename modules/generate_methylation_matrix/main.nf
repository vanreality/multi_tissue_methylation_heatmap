process GENERATE_METHYLATION_MATRIX {
    tag "$meta.id"
    
    input:
    val meta
    path bedgraph_list
    path bed
    path script

    output:
    tuple val(meta), path("*.tsv"), emit: methylation_matrix

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 ${script} \\
        --bedgraph_list ${bedgraph_list} \\
        --bed ${bed} \\
        --output ${prefix}.tsv
    """
}
