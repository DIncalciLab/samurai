def VERSION = "0.1"

process CORRECT_LOGR_ICHORCNA {

    tag "Correcting Log2 for GISTIC Analysis"
    label "process_low"
// TO DO: Create a container in the repository to be pulled
    container "quay.io/dincalcilab/tidyverse:1.0.0-673997e"

    input:
        path(seg_file)
        path(ploidy_summary)


    output:
        path("*_logR_corrected_gistic.seg"),                            emit: gistic_file
        path("versions.yml"),                                           emit: versions

    script:

    """
    correct_logR_ichorcna.R \\
        --seg ${seg_file} \\
        --ploidy ${ploidy_summary}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        correct_logR_ichorcna: ${VERSION}
    END_VERSIONS
    """

}
