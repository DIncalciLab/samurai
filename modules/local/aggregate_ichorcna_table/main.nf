def VERSION = '0.0.1'

process AGGREGATE_ICHORCNA_TABLE {

    tag "Assemble IchorCNA outputs"
    container "quay.io/einar_rainhart/pandas-pandera:1.5.3"
    label 'process_low'

    input:
    file(all_params)


    output:
        path("ichorcna_summary_mqc.txt"),       emit: ichorcna_summary
        //path("ichorcna_data_mqc.tsv"),        emit: ichorcna_summary
        path("ichorCNA_summary.xlsx")
        path "versions.yml",                  emit: versions

    script:

    """
    assemble_outputs.py ${all_params}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        assemble_outputs: ${VERSION}
    END_VERSIONS

    """
}
