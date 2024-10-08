def VERSION = "1.3.0"

process QDNASEQ {

    tag "${meta.id}"
    label "process_medium"
    stageInMode "link"

    container "quay.io/dincalcilab/qdnaseq:1.30.0-a28ebc1"

    input:
        tuple val(meta), path(bamfiles), path(bamindex)
        val(binsize)
        val(genome)

    output:
        path("*_bins.bed"),                       emit: bins, optional: true
        path("*.calls.seg"),                      emit: called_segments
        path("*_.seg"),                           emit: segments
        path("*_filt.seg"),                       emit: filtered_segments
        path("*.rds"),                            emit: processed_data
        path("*_bin_plot.pdf"),                   emit: bin_plot
        path("*_segment_plot.pdf"),               emit: segment_plot, optional: true
        path("*_summary.txt"),                    emit: summary_table
        path("versions.yml"),                     emit: versions

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    run_qdnaseq.R \\
        --cpus "${task.cpus}" \\
        --project-name "${prefix}" \\
        --paired-end ${params.qdnaseq_paired_ends} \\
        --bin-size ${binsize} \\
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
