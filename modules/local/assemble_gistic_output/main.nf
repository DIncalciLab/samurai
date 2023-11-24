VERSION = '0.1'

process ASSEMBLE_GISTIC_OUTPUT {

    tag "Rearrange GISTIC2 output"
    label "process_single"
    //TO DO: Create container in quay.io
    container "/home/sarap/cache_singularity/rearrange_gistic_out.img"

    input:
        path gistic_file


    output:
        path("gistic_lesions.txt"),             emit: gistic_lesions
        // path("gistic_genes.txt"),               emit: gistic_genes NEED TO FIX THIS
        path("gistic_log2R.txt"),               emit: gistic_log_r
        path("gistic_cn_states.txt"),           emit: gistic_cn_states
        path("gistic_lesions_mqc.tsv"),         emit: gistic_lesions_mqc
        path("versions.yml"),                   emit: versions

    script:
    
    """
    assemble_gistic_output.py -gistic_file ${gistic_file}

    cat <<-MULTIQC_HEADER > gistic_lesions_mqc.tsv
    #id: 'gistic_lesions'
    #section_name: 'GISTIC2'
    #description: 'This table shows significant lesions identified by GISTIC2 analysis.'
    #format: 'tsv'
    #plot_type: 'table'
    #headers:
    #  chromosome:
    #    title: 'Chromosome'
    #    placement: 1
    #  cytoband:
    #    title: 'Cytoband'
    #    placement: 2
    #  start:
    #    title: 'Start'
    #    placement: 3
    #    format: '{:,.0f}'
    #  end:
    #    title: 'End'
    #    placement: 4
    #    format: '{:,.0f}'
    #  length:
    #    title: 'Length'
    #    placement: 5
    #    format: '{:,.0f}'
    #  kind:
    #    title: 'CN Kind'
    #    placement: 6
    #  frequency:
    #    title: 'Lesion Frequency'
    #    placement: 7
    #    format: '{:,.2f}'
    #  counts:
    #    title: 'counts'
    #    placement: 8
    #  lesion_type:
    #    title: 'Type of Lesion'
    #    placement: 9
    MULTIQC_HEADER

        cat gistic_lesions.txt >> gistic_lesions_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Assemble Gistic Output: ${VERSION}
    END_VERSIONS
    """


}
