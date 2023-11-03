def VERSION = "0.1"

process QUANTIFY_CIN_SIGNATURES {

    tag "Extracting Signatures"
    label "process_single"

    container "/home/sarap/cache_singularity/CINSigQuant.sif"

    input:
        path(segmentation_file)


    output:
        path("*_plot_by_component.pdf"),          emit: signatures_heatmap
        path("*_plot_activities.pdf"),            emit: sig_activity_plot
        path("*_signatures.rds"),                 emit: signature_rds
        path("*_platinum_prediction.rds"),        emit: platinum_prediction
        path("*_activity.txt"),                   emit: sig_activity_file
        path("versions.yml"),                     emit: versions

    script:

    def args = task.ext.args ?: ''

    """
    compute_signatures.R \\
        --seg_file ${segmentation_file} \\
        --cpus "${task.cpus}" \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CINSignatureQuantification: ${VERSION}
    END_VERSIONS
    """

}
