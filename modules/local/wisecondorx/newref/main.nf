def VERSION = "1.2.4"

process WISECONDORX_NEWREF {

    tag "reference"
    label "process_high"
    container "quay.io/biocontainers/wisecondorx:1.2.5--pyh5e36f6f_0"

    input:
        path npzfiles
    output:
        path("*.npz"), emit: npz_reference
        path "versions.yml", emit: versions

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "reference"

    """
    WisecondorX newref \\
        --cpus ${task.cpus} \\
        $args \\
        ${npzfiles} \\
        ${prefix}.npz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wisecondorX: ${VERSION}
    END_VERSIONS
    """

}
