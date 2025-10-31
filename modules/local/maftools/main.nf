process MAFTOOLS {

    tag "Maftools Plot"
    label "process_low"

    container "quay.io/dincalcilab/maftools:2.17.0-ce17bf9"

    input:
    path all_lesions
    path amp_genes
    path del_genes
    path gistic_scores

    output:
    path ("maftools_summary_mqc.png"), emit: chrom_plot
    path ("maftools_bubble.png"), emit: bubble_plot
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    def VERSION = "0.1"

    """
    run_maftools.R   \\
        --all_lesions ${all_lesions} \\
        --amp_genes ${amp_genes} \\
        --del_genes ${del_genes} \\
        --gistic_scores ${gistic_scores} \\
        ${args}



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MAFTOOLS: ${VERSION}
    END_VERSIONS
    """
}
