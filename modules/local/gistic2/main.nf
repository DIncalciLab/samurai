def VERSION="0.1"
process GISTIC2 {

    tag "Run Gistic2"
    label "process_low"

    container "quay.io/dincalcilab/gistic:2.0.23-35a282"

    input:
        path(seg_file)


    output:
        path("all_lesions.conf_*"),               emit: all_lesions
        path("amp_genes.conf*"),                  emit: amplified_genes
        path("del_genes.conf*"),                  emit: deleted_genes
        path("scores.gistic"),                    emit: gistic_score
        path("raw_copy_number.pdf"),              emit: raw_copy_number
        path("amp_qplot.pdf"),                    emit: amp_score_qplot
        path("del_qplot.pdf"),                    emit: del_score_qplot
        path("gistic_results", type: "dir"),      emit: gistic_results_dir
        path("versions.yml"),                     emit: versions

    script:

    def args = task.ext.args ?: ''
    // Inside the container
    def ref_gene_file = params.genome == 'hg38' ? "-refgene '/opt/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat'" : "-refgene '/opt/refgenefiles/hg19.UCSC.add_miR.140312.refgene.mat'"

    """
    gistic2 \\
        -seg '${seg_file}' \\
        ${ref_gene_file} \\
        -ta ${params.gistic_t_amp} \\
        -td ${params.gistic_t_del} \\
        -rx ${params.gistic_remove_x} \\
        -conf ${params.gistic_conf} \\
        -b ./

    mkdir -p gistic_results
    cp all_lesions.conf* {amp,del}_genes.conf* scores.gistic *.pdf  gistic_results/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GISTIC2: ${VERSION}
    END_VERSIONS
    """
}
