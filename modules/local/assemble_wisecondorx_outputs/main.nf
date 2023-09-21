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
        path "wisecondorx_summary_mqc.txt"   
        path "wisecondorx_summary_mqc.tsv"   , emit: wisecondorx_summary
        path "all_aberrations.txt"           , emit: aberration_stats
        path "versions.yml"                  , emit: versions

    script:

    """

    assemble_statistics.py ${statistics}
    assemble_aberrations.py ${aberrations}

    cat <<-MULTIQC_HEADER > wisecondorx_summary_mqc.tsv
    #id: 'wisecondorx_data_summary'
    #section_name: 'WisecondorX'
    #description: 'This table shows the Copy Number Score (CPA) and Median Segment Variance (MSV) computed by WisecondorX.'
    #format: 'tsv'
    #plot_type: 'table'
    #headers:
    #  sample:
    #    title: 'Sample ID'
    #    placement: 1
    #  reads:
    #    title: 'Number of Reads'
    #    placement: 2
    #  cpa:
    #    title: 'Copy number profile abnormality Score'
    #    format: '{:,.4f}'
    #    placement: 3
    #    scale: False
    #  msv:
    #    title: 'Median segment variance'
    #    format: '{:,.4f}'
    #    placement: 4
    #    scale: False
    MULTIQC_HEADER

    # Reorder samples
    cat wisecondorx_summary_mqc.txt >> wisecondorx_summary_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        assemble_statistics: ${VERSION}
        assemble_aberrations: ${VERSION}
    END_VERSIONS

    """


}