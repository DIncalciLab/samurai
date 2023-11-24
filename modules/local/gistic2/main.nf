def VERSION="0.1"
process GISTIC2 {

    tag "Run Gistic2"
    label "process_low"

    // TO DO: Create container on quay.io
    container '/mnt/svgs/cache_singularity/gistic.img'

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
        path("versions.yml"),                     emit: versions

    script:

    def args = task.ext.args ?: ''
    def ref_gene_file = params.genome == 'hg38' ? "-refgene '/opt/GISTIC/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat'" : "-refgene '/opt/GISTIC/refgenefiles/hg19.UCSC.add_miR.140312.refgene.mat'"
    
    """     
    /usr/local/bin/gistic2 \\
        -seg '${seg_file}' \\
        ${ref_gene_file} \\
        -ta ${params.gistic_t_amp} \\
        -td ${params.gistic_t_del} \\
        -rx ${params.gistic_remove_x} \\
        -conf ${params.gistic_conf} \\
        -b ./


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GISTIC2: ${VERSION}
    END_VERSIONS
    """
}
