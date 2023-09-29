def VERSION = '0.0.1'

process CREATE_ASCATSC_SUMMARY{

    tag "create ASCAT.sc summary for multiqc"
    label "process_single"

    input:
        path "ascatsc_summary.txt"
        
    output:
        path "ascatsc_summary_mqc.tsv"       , emit: ascatsc_summary
        path "versions.yml"                  , emit: versions

    script:

        """
        cat <<-MULTIQC_HEADER > ascatsc_summary_mqc.tsv
        #id: 'wisecondorx_data_summary'
        #section_name: 'ASCATsc'
        #description: 'This table shows purity and ploidy values computed by ASCAT.sc.'
        #format: 'tsv'
        #plot_type: 'table'
        #headers:
        #  samplename:
        #    title: 'Sample ID'
        #    placement: 1
        #  purity:
        #    title: 'Purity'
        #    format: '{:,.4f}'
        #    placement: 2
        #    scale: False
        #  ploidy:
        #    title: 'Ploidy'
        #    format: '{:,.4f}'
        #    placement: 3
        #    scale: False
        MULTIQC_HEADER

        # Reorder samples
        cat ascatsc_summary.txt >> ascatsc_summary_mqc.tsv
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            create_ascatsc_summary: ${VERSION}
        END_VERSIONS
        """
}