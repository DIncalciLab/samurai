def VERSION = "1.2.4"

process WISECONDORX_CONVERT {

    tag "${meta.id}"
    label "process_low"
    container "quay.io/biocontainers/wisecondorx:1.2.5--pyh5e36f6f_0"

    input:
        tuple val(meta), path(bamfile), path(bamindex)
    output:
        tuple val(meta), path("*.npz"), emit: npz
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    WisecondorX convert \\
        ${args} \\
        ${bamfile} \\
        ${prefix}.npz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wisecondorX: ${VERSION}
    END_VERSIONS
    """

}