def VERSION="0.1"
process RUN_GISTIC2 {

    tag "Run Gistic2"
    label "process_low"
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
        -ta ${params.gistic_ta} \\
        -td ${params.gistic_td} \\
        -rx ${params.gistic_rx} \\
        -conf ${params.gistic_conf} \\
        -b ./

    cat <<-MULTIQC_HEADER > mqc_amp_qplot.pdf
    id: 'amplification_qplot'
    section_name: 'GISTIC Amplification'
    description: 'This plot shows amplified genes identified by GISTIC2.'
    format: 'png'
    MULTIQC_HEADER

    cat amp_qplot.pdf >> mqc_amp_qplot.pdf
    

    cat <<-MULTIQC_HEADER > mqc_del_qplot.pdf
    id: 'deletion_qplot'
    section_name: 'GISTIC Deletion'
    description: 'This plot shows deleted genes identified by GISTIC2.'
    format: 'png'
    MULTIQC_HEADER

    cat del_qplot.pdf >> mqc_del_qplot.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GISTIC2: ${VERSION}
    END_VERSIONS
    """
}