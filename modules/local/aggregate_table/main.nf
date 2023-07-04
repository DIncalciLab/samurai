def VERSION = '0.0.1'

process AGGREGATE_TABLE {

    container "quay.io/einar_rainhart/pandas-pandera:1.5.3"

    input:
    file(all_params)


    output:
        path("ploidy_summary.txt"), emit: ploidy_summary
        path("ichorCNA_summary.xlsx")
        path "versions.yml", emit: versions

    script:

    """
    assemble_outputs.py ${all_params}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        assemble_outputs: ${VERSION}
    END_VERSIONS

    """
}
