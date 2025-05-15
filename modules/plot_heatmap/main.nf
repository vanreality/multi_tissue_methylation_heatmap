process PLOT_HEATMAP {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(methylation_matrix)
    path script

    output:
    tuple val(meta), path("*.pdf"), emit: heatmap

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ""
    """
    python3 ${script} \\
        --input ${methylation_matrix} \\
        --output-plot ${prefix}.pdf \\
        ${args}
    """
}
