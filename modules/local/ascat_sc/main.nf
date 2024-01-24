def VERSION = "0.1"

process ASCAT_SC {

    tag "${meta.id}"
    label "process_medium"

    container "quay.io/dincalcilab/ascat_sc:1.0.0"

    input:
        tuple val(meta), path(bamfiles), path(bamindex)
        val(binsize)
        val(genome)


    output:
        path("profiles_${meta.id}.pdf"),                       emit: profiles_plot
        path("profiles_${meta.id}_refitted.pdf"),              emit: profiles_refitted
        //path("*_gistic.seg"),                     emit: gistic_file
        path("*_summary.txt"),                    emit: summary_table
        path("*_segments.seg"),                   emit: segments
        path("*.rds"),                            emit: rds
        path("*_df_signatures.seg"),              emit: sig_file
        path("*.rds"),                            emit: ascat_rds
        path("versions.yml"),                     emit: versions

    script:

    def args = task.ext.args ?: ''
    def gender = meta.gender ? "--sex ${meta.gender}": ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bins_bp = binsize ? "--binsize ${binsize}000": "--binsize 30000"

    """
    run_ascatsc.R \\
        --tumor_bams ${bamfiles} \\
        --cpus "${task.cpus}" \\
        --projectname "${prefix}" \\
        --genome "${genome}" \\
        ${bins_bp} \\
        ${gender} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ASCAT.sc: ${VERSION}
    END_VERSIONS
    """

}
