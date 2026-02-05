process FORMAT_ICHORCNA_SEG {
    tag "$meta.id"
    container "quay.io/einar_rainhart/pandas-pandera:1.5.3"
    label 'process_low'

    input:
    tuple val(meta), path(seg)

    output:
    path "${meta.id}_formatted.seg", emit: seg

    script:
    def VERSION = '0.0.1'

    """
    reformat_ichorcna_seg.py -s $seg --o ${meta.id}_formatted.seg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        format_ichorcna: ${VERSION}
    END_VERSIONS
    """
}
