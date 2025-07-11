process PLOT_ICHORCNA {

    tag "Plotting ichorCNA results for $meta.id"
    label 'process_low'

    container "community.wave.seqera.io/library/procps-ng_r-argparser_r-dplyr_r-ggplot2_pruned:10da72fa04bcba1a"

    input:
        tuple val(meta), path(cna_seg)
        tuple val(meta), path(bins)
        file(ichorcna_params)

    output:
        tuple val(meta), path("*.copy_number.png")        , emit: plot_png
        tuple val(meta), path("*.copy_number.svg")        , emit: plot_svg
        path "versions.yml"                               , emit: versions

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def VERSION = '0.1'

        """
        plot_ichorcna_output.R \\
            --id ${prefix} \\
            --seg_file ${cna_seg} \\
            --binfile ${bins} \\
            --summary ${ichorcna_params} \\
            --outdir .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            plot_ichorcna: $VERSION
        END_VERSIONS
        """
}
