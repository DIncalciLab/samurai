def VERSION = "1.2.4"

process WISECONDORX_PREDICT {

    tag "${meta.id}"
    label "process_single"
    container "quay.io/biocontainers/wisecondorx:1.2.5--pyh5e36f6f_0"

    input:
        tuple val(meta), path(npzfile)
        path(blacklist)
        path(reference_data)
    output:
        tuple val(meta), path("*_bins.bed"),            emit: bins
        tuple val(meta), path("*_aberrations.bed"),     emit: calls
        tuple val(meta), path("*_segments.bed"),        emit: segments
        tuple val(meta), path("*_statistics.txt"),      emit: statistics
        tuple val(meta), path("*.plots"),               emit: plots
        tuple val(meta), path("*_genome_plot.png"),     emit: genome_plot
        path "versions.yml",                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def blacklist_bed = blacklist ? "--blacklist ${blacklist}" : ''
    def gender = meta.gender ? "--gender ${meta.gender}": ''

    """
    WisecondorX predict \\
        ${npzfile} \\
        ${reference_data} \\
        ${prefix} \\
        $args \\
        ${blacklist_bed} \\
        ${gender} \\
        --plot \\
        --bed \\
        --add-plot-title

    cp ${prefix}.plots/genome_wide.png ${prefix}_genome_plot.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wisecondorX: ${VERSION}
    END_VERSIONS
    """

}