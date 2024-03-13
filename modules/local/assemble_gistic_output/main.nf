VERSION = '0.1'

process ASSEMBLE_GISTIC_OUTPUT {

    tag "Rearrange GISTIC2 output"
    label "process_single"
    container "quay.io/dincalcilab/gistic-cli:0.4.2-59c51ab"

    input:
        path gistic_folder

    output:
        path("gistic_summary_mqc.txt"),             emit: gistic_lesions
        path("gistic.gene_table.txt"),              emit: gistic_genes
        path("gistic_log2R.txt"),                   emit: gistic_log_r
        path("gistic_cn_states.txt"),               emit: gistic_cn_states
        path("versions.yml"),                       emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''

    """
    gistic lesion_table ${args} ${gistic_folder} gistic_summary_mqc
    mv gistic_summary_mqc.focal_lesions.txt gistic_summary_mqc.txt

    gistic gene_table ${args} ${gistic_folder} gistic
    gistic copy_number_state_table ${args} ${gistic_folder} gistic_cn_states.txt
    gistic log2ratio_table ${args} ${gistic_folder} gistic_log2R.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Assemble Gistic Output: ${VERSION}
    END_VERSIONS
    """


}
