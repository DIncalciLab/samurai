process COMPUTE_CINMETRICS{
    tag "Compute CINmetrics"
    label 'process_single'

    container 'ghcr.io/dincalcilab/cinmetrics:0.1.0'

    input:
        file(seg_file)

    output:
        path("*.tsv"),              emit: cinmetrics_summary
        path("versions.yml"),       emit: versions

    when:
    task.ext.when == null || task.ext.when

   script:
    def args = task.ext.args ?: ''
    def VERSION = "0.1"
    """
    compute_cinmetrics.R \\
        --seg_file $seg_file \\
        --output cinmetrics_summary_mqc.tsv \\
$args

    cat << END_VERSIONS > versions.yml
    "${task.process}":
        CINmetrics: ${VERSION}
END_VERSIONS
    """
}