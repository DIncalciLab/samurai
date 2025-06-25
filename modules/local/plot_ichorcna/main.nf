process PLOT_ICHORCNA {

    tag "Plotting ichorCNA results for $meta.id"
    label 'process_low'

    container "/mnt/local/sarap/images_def/plot_ichorcna/plot_ichorcna.sif"

    input:
        tuple val(meta), path(cna_seg)
        tuple val(meta), path(bins)
        file(ichorcna_summary)

    output:
        tuple val(meta), path("*.copy_number.png")        , emit: corrected_cn_plot
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
            --summary ${ichorcna_summary} \\
            --outdir .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            plot_ichorcna: $VERSION
        END_VERSIONS
        """
}
