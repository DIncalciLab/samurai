def VERSION = "0.1"

process CORRECT_LOGR_ICHORCNA {

    tag "${meta.id}"
    label "process_single"

    container "/home/sarap/cache_singularity/dplyr_readr.sif"

    input:
        tuple val(meta), path(seg_file)
        path(ploidy_summary)


    output:
        tuple val(meta), path("*.seg"),           emit: gistic_file
        path("versions.yml"),                     emit: versions
        
    script:

    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    correct_logR_ichorcna.R \\
        --seg ${seg_file} \\
        --ploidy ${ploidy_summary} \\
        --project ${prefix} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        correct_logR_ichorcna: ${VERSION}
    END_VERSIONS
    """

}
