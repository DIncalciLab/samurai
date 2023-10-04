def VERSION = "0.1"

process ASCAT_SC {

    tag "${meta.id}"
    label "process_medium"

    container "/home/sarap/cache_singularity/ASCAT-sc.sif"

    input:
        tuple val(meta), path(bamfiles), path(bamindex)


    output:
        path("profiles_*"),                       emit: profiles_plot
        path("*_gistic.seg"),                     emit: gistic_file
        path("*_summary.txt"),                    emit: summary_table
        path("*_segments.seg"),                   emit: segments
        path("versions.yml"),                     emit: versions

    script:

    def args = task.ext.args ?: ''
    def gender = meta.gender ? "--sex ${meta.gender}": ''
         
    """
    run_ascatsc.R \\
        --tumour_bams ${bamfiles} \\
        --cpus "${task.cpus}" \\
        --projectname "${meta.id}" \\
        ${gender} \\
        $args 
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ASCAT.sc: ${VERSION}
    END_VERSIONS
    """

}
