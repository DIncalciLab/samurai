process PLOT_WISECONDORX_CNV {
    tag "Plotting WisecondorX results for $meta.id"
    label 'process_low'
    container "community.wave.seqera.io/library/procps-ng_r-argparser_r-dplyr_r-ggplot2_pruned:10da72fa04bcba1a"

    input:
        tuple val(meta), path(seg_file)
        tuple val(meta2), path(bins)

    output:
        tuple val(meta), path("*.copy_number.png"), emit: plot_png
        tuple val(meta), path("*.copy_number.svg"), emit: plot_svg
        path "versions.yml"                       , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ''
        def prefix  = task.ext.prefix ?: "${meta.id}"
        def VERSION = '0.1'
        """
        plot_wisecondorx_cnv.R \\
            --id ${prefix} \\
            --seg_file ${seg_file} \\
            --binfile ${bins} \\
            --outdir . \\
            ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            plot_wisecondorx_cnv: $VERSION
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.copy_number.png
        touch ${prefix}.copy_number.svg
        touch versions.yml
        """
}
