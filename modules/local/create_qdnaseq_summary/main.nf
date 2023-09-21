def VERSION = '0.0.1'

process CREATE_QDNASEQ_SUMMARY{

    tag "create qdnaseq summary for multiqc"
    label "process_single"

    input:
        path "qdnaseq_summary.txt"
        
    output:
        path "qdnaseq_summary_mqc.tsv"       , emit: qdnaseq_summary
        path "versions.yml"                  , emit: versions

    script:

        """
        cat <<-MULTIQC_HEADER > qdnaseq_summary_mqc.tsv
        #id: 'qdnaseq_data_summary'
        #section_name: 'QDNAseq'
        #description: 'This table shows a summary table of QDNAseq analysis.'
        #format: 'tsv'
        #plot_type: 'table'
        #headers:
        #  sample:
        #    title: 'Sample ID'
        #    placement: 1
        #  nr_reads:
        #    title: 'Number of Reads'
        #    placement: 2
        #  expected_var:
        #    title: 'Expected Variance
        #    format: '{:,.4f}
        #    placement: 3
        #    scale: False
        #  expected_sd:
        #    title: 'Expected Standard Deviation'
        #    format: '{:,.4f}'
        #    placement: 4
        #    scale: False
        #  binsize:
        #    title: 'Size of bins'
        #    format: '{:,.0f}'
        #    placement: 5
        MULTIQC_HEADER

        cat qdnaseq_summary.txt >> qdnaseq_summary_mqc.tsv
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            create_qdna_summary: ${VERSION}
        END_VERSIONS
        """
}