process ASSEMBLE_GISTIC_OUTPUT {

    tag "Rearrange GISTIC2 output"
    label "process_single"
    container "ghcr.io/dincalcilab/gistic-cli:0.5.4-ade3d88"

    input:
        path gistic_folder

    output:
        path("gistic_summary_mqc.txt"),                emit: gistic_lesions
        path("gistic.gene_table.txt"),                 emit: gistic_genes
        path("gistic_log2R.txt"),                      emit: gistic_log_r
        path("gistic_cn_states.txt"),                  emit: gistic_cn_states
        path("gistic_broad_lesions_mqc.txt"),          emit: gistic_broad_lesions, optional: true
        path("versions.yml"),                          emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''

    """
    gistic lesion_table -t all ${args} ${gistic_folder} gistic_summary_mqc
    mv gistic_summary_mqc.focal_lesions.txt gistic_summary_mqc.txt

    # Only do the rename if the file was produced, i.e. when broad analysis runs
    if [ -f gistic_summary_mqc.broad_lesions.txt ]; then
        mv gistic_summary_mqc.broad_lesions.txt gistic_broad_lesions_mqc.txt
    fi

    gistic gene_table ${args} ${gistic_folder} gistic
    gistic copy_number_state_table ${args} ${gistic_folder} gistic_cn_states.txt
    gistic log2ratio_table ${args} ${gistic_folder} gistic_log2R.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gistic-cli: \$(gistic --version | cut -d ' ' -f2)
    END_VERSIONS
    """

}
