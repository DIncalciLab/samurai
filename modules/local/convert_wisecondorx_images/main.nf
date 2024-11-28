process CONVERT_WISECONDORX_IMAGES {

    tag "genome_plot"
    container "community.wave.seqera.io/library/pip_img2pdf:75fb1abd67dc4e95"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_img2pdf:43f78f9b1d6c437d' :
        'community.wave.seqera.io/library/pip_img2pdf:75fb1abd67dc4e95' }"

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
