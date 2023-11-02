process CONCATENATE_SEG {

    tag "Concatenate Segment files"
    label 'cpu_low'

    container "quay.io/dincalcilab/qdnaseq:1.30.0"

    input:
        path(seg_files)
        val output_name
    output:
        path("*.seg"), emit: all_segments

    script:

    """
    concatenate_segments.R \\
        --segfiles ${seg_files} \\
        --output ${output_name}
    """

}
