process CONCATENATE_PDF {

    tag "Concatenate_PDF"
    label 'cpu_low'

    container 'docker.io/t0shy/qpdf-docker:11.3.0'

    input:
        path pdf_files
        
    output:
        path("genome_plot.pdf"), emit: genome_plot
        path "versions.yml", emit: versions

    script:

    """
    qpdf --empty --pages ${pdf_files} -- genome_plot.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qpdf: \$(qpdf --version | sed -n '1p' | cut -d' ' -f3)
    END_VERSIONS
    """

}