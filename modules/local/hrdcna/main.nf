def VERSION = "0.1"

process HRDCNA {

    tag "Computing HRDCNA score"
    label "process_single"

    container "quay.io/dincalcilab/hrdcna:0.0.0.9000_deb8861"

    input:
        file(segmentation_file)

    output:
        path("hrdcna_summary_mqc.tsv"),             emit: hrdcna_summary
        path("hrdcna_scores.rds"),                  emit: hrdcna_scores_rds
        path("hrdcna_features_activity.tsv"),       emit: features_activity
        path("hrdcna_features_activity.rds"),       emit: features_activity_rds
        path("versions.yml"),                       emit: versions

    when:
    task.ext.when == null || task.ext.when
    
    script:

    def args = task.ext.args ?: ''

    """
    run_hrdcna.R \\
        --seg_file '${segmentation_file}' \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        HRDCNA: ${VERSION}
    END_VERSIONS
    """




}       