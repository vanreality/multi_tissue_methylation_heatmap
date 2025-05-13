process CALCULATE_DMR_METHYLATION {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(bedgraph)
    path script

    output:
    tuple val(meta), path("*_DMR.bedGraph"), emit: bedgraph

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 ${script} \\
        --bedgraph ${bedgraph} \\
        --output ${prefix}_DMR.bedGraph
    """
}
