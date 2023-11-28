def VERSION = '0.0.1'

process AGGREGATE_ICHORCNA_TABLE {

    tag "Assemble IchorCNA outputs"
    container "quay.io/einar_rainhart/pandas-pandera:1.5.3"
    label 'process_low'

    input:
    file(all_params)


    output:
        path("ploidy_summary_mqc.txt"),       emit: ploidy_summary
        path("ichorcna_data_mqc.tsv"),        emit: ichorcna_summary
        path("ichorCNA_summary.xlsx")
        path "versions.yml",                  emit: versions

    script:

    """
    assemble_outputs.py ${all_params}

    cat <<-MULTIQC_HEADER > ichorcna_data_mqc.tsv
    #id: 'ichorcna_data_summary'
    #section_name: 'ichorCNA'
    #description: 'This table shows the tumor fraction estimates made by ichorCNA.'
    #format: 'tsv'
    #plot_type: 'table'
    #col1_header: 'Sample ID'
    #headers:
    #  samplename:
    #    title: 'Sample ID'
    #    placement: 1
    #  Ploidy:
    #    title: 'Ploidy'
    #    format: '{:,.0f}'
    #    placement: 2
    #    scale: False
    #  'Tumor Fraction':
    #    title: 'Tumor fraction'
    #    format: '{:,.4f%}'
    #    placement: 3
    #    scale: False
    #  'GC-Map correction MAD':
    #    title: 'MAD'
    #    format: '{:,.4f}'
    #    scale: false
    #    placement: 4
    #    max: 0.5
    #    cond_formatting_rules:
    #      pass:
    #        - lt: 0.15
    #      warn:
    #        - gt: 0.15
    #      fail:
    #        - gt: 0.25
    MULTIQC_HEADER

    # Reorder samples
    cat ploidy_summary_mqc.txt >> ichorcna_data_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        assemble_outputs: ${VERSION}
    END_VERSIONS

    """
}
