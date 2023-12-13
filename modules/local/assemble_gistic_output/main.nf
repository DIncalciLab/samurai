VERSION = '0.1'

process ASSEMBLE_GISTIC_OUTPUT {

    tag "Rearrange GISTIC2 output"
    label "process_single"
    //TO DO: Create container in quay.io
    container "/home/sarap/cache_singularity/rearrange_gistic_out.img"

    input:
        path gistic_file

    output:
        path("gistic_lesions_mqc.txt"),             emit: gistic_lesions
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
