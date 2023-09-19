def VERSION = '0.0.1'

process ASSEMBLE_WISECONDORX_OUTPUTS {

    tag "output tables"
    container "quay.io/einar_rainhart/pandas-pandera:1.4.3"
    label "process_single"

    input:
        path statistics
        path aberrations
    output:
        path "summary.xlsx"                  , emit: summary
        path "*summary.txt"                  , emit: wisecondorx_summary
        path "all_aberrations.txt"           , emit: aberration_stats
        path "versions.yml"                  , emit: versions

    script:

    """

    assemble_statistics.py ${statistics}
    assemble_aberrations.py ${aberrations}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        assemble_statistics: ${VERSION}
        assemble_aberrations: ${VERSION}
    END_VERSIONS

    """


}