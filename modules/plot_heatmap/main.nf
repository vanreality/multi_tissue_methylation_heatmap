process PLOT_HEATMAP {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(methylation_matrix)
    path script

    output:
    tuple val(meta), path("*.pdf"), emit: heatmap

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 ${script} \\
        --methylation_matrix ${methylation_matrix} \\
        --output ${prefix}
    """
}
