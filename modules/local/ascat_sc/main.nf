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
        path("profiles_*"),                       emit: profiles_plot
        path("*_gistic.seg"),                     emit: gistic_file
        path("*_summary.txt"),                    emit: summary_table
        path("*_segments.seg"),                   emit: segments
        path("*.rds"),                            emit: rds
        path("versions.yml"),                     emit: versions

    script:

    def args = task.ext.args ?: ''
    def gender = meta.gender ? "--sex ${meta.gender}": ''
    def prefix = task.ext.prefix ?: "${meta.id}"


    """
    run_ascatsc.R \\
        --tumour_bams ${bamfiles} \\
        --cpus "${task.cpus}" \\
        --projectname "${prefix}" \\
        --binsize "${binsize}" \\
        --genome "${genome}" \\
        ${gender} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ASCAT.sc: ${VERSION}
    END_VERSIONS
    """

}
