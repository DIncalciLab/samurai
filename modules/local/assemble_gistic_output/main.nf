VERSION = '0.1'

process ASSEMBLE_GISTIC_OUTPUT {

    tag "Rearrange GISTIC2 output"
    label "process_single"
    container "quay.io/dincalcilab/gistic-cli:0.4.1-3e0259a"

    input:
        path gistic_file

    output:
        path("gistic_summary_mqc.txt"),             emit: gistic_lesions
        // path("gistic_genes.txt"),               emit: gistic_genes NEED TO FIX THIS
        path("gistic_log2R.txt"),               emit: gistic_log_r
        path("gistic_cn_states.txt"),           emit: gistic_cn_states
        path("versions.yml"),                   emit: versions

    script:

    """
    assemble_gistic_output.py -gistic_file ${gistic_file}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Assemble Gistic Output: ${VERSION}
    END_VERSIONS
    """


}
