def VERSION = "0.1"

process CIN_SIGNATURE_QUANTIFICATION {

    tag "Extracting Signatures"
    label "process_single"

    container "quay.io/dincalcilab/cinsignaturequantification:1.1.2-af6117b"

    input:
        file(segmentation_file)


    output:
        path("*_plot_by_component.png"),          emit: signatures_heatmap
        path("signatures_summary_mqc.png"),       emit: sig_activity_plot
        path("*_signatures.rds"),                 emit: signature_rds
        path("*_activity.txt"),                   emit: sig_activity_file
        path("versions.yml"),                     emit: versions

    script:

    def args = task.ext.args ?: ''

    """
    compute_signatures.R \\
        --seg_file '${segmentation_file}' \\
        --cpus '${task.cpus}' \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CINSignatureQuantification: ${VERSION}
    END_VERSIONS
    """

}
