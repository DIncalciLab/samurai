process CONVERT_GISTIC_SEG {

    tag "${meta.id}"
    container "quay.io/dincalcilab/pandas-pybedtools:1.4.4-ce17bf9"
    label "process_single"

    input:
        tuple val(meta), path(segments)
        tuple val(meta), path(bins)
        tuple val(meta), path(aberrations)
    output:
        tuple val(meta), path("*.seg"), emit: segfile
        tuple val(meta), path("*_gistic.seg"), emit: gistic_file
        path "versions.yml",           emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:

    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    convert_segment.py \\
        --id ${prefix} \\
        ${segments} \\
        ${bins} \\
        ${aberrations} \\
        ${prefix}.seg \\
        ${prefix}_gistic.seg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        convert_segment: \$(convert_segment.py --version | cut -d' ' -f2)
    END_VERSIONS

    """

}
