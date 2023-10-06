process CONVERT_WISECONDORX_IMAGES {

    tag "genome_plot"
    container "quay.io/einar_rainhart/img2pdf:latest"
    label "process_single"

    input:
        path images
    output:
        path "genome_plots.pdf", emit: genome_plot
        path "versions.yml",     emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """

        img2pdf ${images} \\
            --title "WisecondorX results" \\
           -o genome_plots.pdf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            img2pdf: \$(img2pdf --version | cut -d' ' -f2)
        END_VERSIONS

        """
}
