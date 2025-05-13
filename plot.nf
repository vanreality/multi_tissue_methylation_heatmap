include { BEDTOOLS_EXTRACT_SITES } from './modules/bedtools_extract_sites/main.nf'
include { CALCULATE_DMR_METHYLATION } from './modules/calculate_dmr_methylation/main.nf'
include { GENERATE_METHYLATION_MATRIX } from './modules/generate_methylation_matrix/main.nf'

workflow {
    // Read input CSV file and create a channel with sample ID and bedgraph file path
    ch_input = Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row -> 
            def meta = [id: row.sample]
            [ meta, file(row.bedgraph_file_path) ]
        }

    // Pass the channel to BEDTOOLS_EXTRACT_SITES process
    BEDTOOLS_EXTRACT_SITES(ch_input, params.input_bed)

    // Pass the channel to CALCULATE_DMR_METHYLATION process
    CALCULATE_DMR_METHYLATION(
        BEDTOOLS_EXTRACT_SITES.out, 
        file("${workflow.projectDir}/bin/calculate_dmr_methylation.py")
    )

    // Collect all CALCULATE_DMR_METHYLATION outputs and create a CSV file
    ch_bedgraph_list = CALCULATE_DMR_METHYLATION.out
        .collectFile(
            name: 'dmr_methylation_results.csv',
            newLine: true,
            seed: "sample,bedgraph_file_path\n"
        ) { meta, bedgraph ->
            def pubPath = "${params.output_dir}/${meta.id}/${bedgraph.getName()}"
            "${meta.id},${pubPath}"
        }

    // Generate methylation matrix from the collected results
    GENERATE_METHYLATION_MATRIX(
        [id: 'methylation_matrix'], 
        ch_bedgraph_list,
        file("${workflow.projectDir}/bin/generate_methylation_matrix.py")
    )
}

