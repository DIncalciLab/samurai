def VERSION = "0.1"

process REARRANGE_ICHORCNA_OUTPUT {

    tag "Rearranging IchorCNA Output"
    label "process_single"

    container "/mnt/svgs/cache_singularity/dplyr_readr.sif"

    input:
        path(segmentation_file)


    output:
        path("*.seg")                                , emit: sig_file
        path "versions.yml"                          , emit: versions

    script:

    """
    rearrange_ichorcna_output.R \\
        --seg_file ${segmentation_file} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rearrange_ichorcna_output: ${VERSION}
    END_VERSIONS
    """

}