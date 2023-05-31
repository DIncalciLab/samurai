def VERSION = "1.3.0"

process QDNASEQ {

    tag "${meta.id}"
    label "process_high"
    stageInMode "link"

    container "quay.io/dincalcilab/qdnaseq:1.30.0"

    input:
        tuple val(meta), path(bamfiles), path(bamindex)
        val(binsizes)

    output:
        path("*_bins.bed"),                       emit: bins, optional: true
        path("*.calls.seg"),                      emit: called_segments
        path("*_.seg"),                           emit: segments
        path("*_filt.seg"),                       emit: filtered_segments
        path("*.rds"),                            emit: processed_data
        path("*_bin_plot.pdf"),                   emit: bin_plot
        path("*_segment_plot.pdf"),               emit: segment_plot, optional: true
        path("*_summary.txt"),                     emit: summary_table
        path("versions.yml"),                     emit: versions

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def binfile = binsizes ? "--bin-data ${binsizes}" : ''

    """
    run_qdnaseq.R \\
        --cpus "${task.cpus}" \\
        ${binfile} \\
        --project-name "${prefix}" \\
        --paired-end \\
        $args \\
        --source ${bamfiles} \\
        -- \\
        ./


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qDNAseq: ${VERSION}
    END_VERSIONS
    """

}
