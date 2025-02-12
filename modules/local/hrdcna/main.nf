def VERSION = "0.1"

process RUN_HRDCNA {

    tag "Computing HRDCNA score"
    label "process_single"

    container "/mnt/svgs/cache_singularity/hrdcna_0.1.sif"

    input:
        file(segmentation_file)

    output:
        path("hrdcna_scores.tsv"),                  emit: hrdcna_scores
        path("hrdcna_scores.rds"),                  emit: hrdcna_scores_rds
        path("hrdcna_features_activity.tsv"),       emit: features_activity
        path("hrdcna_features_activity.rds"),       emit: features_activity_rds
        path("versions.yml"),                       emit: versions

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