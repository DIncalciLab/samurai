def VERSION = "0.1"

process CORRECT_LOGR_ICHORCNA {

    tag "Correcting Log2 for GISTIC Analysis"
    label "process_low"

    container "/home/sarap/cache_singularity/dplyr_readr.sif"

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
