def VERSION="0.1"
process MAFTOOLS {

    tag "Maftools Plot"
    label "process_low"

    container "quay.io/dincalcilab/maftools:2.17.0-ce17bf9"

    input:
        path(all_lesions)
        path(amp_genes)
        path(del_genes)
        path(gistic_scores)


    output:
        path("plot_gistic_chrom_mqc.png"),        emit: chrom_plot
        path("plot_gistic_bubble.png"),           emit: bubble_plot
        path("versions.yml"),                     emit: versions

    script:

    def args = task.ext.args ?: ''

    """
    run_maftools.R   \\
        --all_lesions ${all_lesions} \\
        --amp_genes ${amp_genes} \\
        --del_genes ${del_genes} \\
        --gistic_scores ${gistic_scores} \\
        $args



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MAFTOOLS: ${VERSION}
    END_VERSIONS
    """
}
