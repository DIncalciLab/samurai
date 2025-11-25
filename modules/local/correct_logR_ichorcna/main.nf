process CORRECT_LOGR_ICHORCNA {

    tag "Correcting Log2 for GISTIC Analysis"
    label "process_low"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/83908efbfba116e554b59ef27409fc78f1f2cc94158c5e01dde680527c3de567/data'
        : 'community.wave.seqera.io/library/polars_procps-ng_typer:d1a53d7945a021e3'}"

    input:
    path seg_file
    path ploidy_summary

    output:
    path ("*_logR_corrected_gistic.seg"), emit: gistic_file
    path ("versions.yml"), emit: versions

    script:
    def VERSION = "0.1"

    """
    correct_logr_ichorcna.py \\
        --seg ${seg_file} \\
        --ploidy ${ploidy_summary}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        correct_logR_ichorcna: ${VERSION}
    END_VERSIONS
    """
}
