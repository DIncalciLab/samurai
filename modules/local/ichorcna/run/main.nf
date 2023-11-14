process RUN_ICHORCNA {

    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::r-ichorcna=0.3.2"
    container "/mnt/svgs/cache_singularity/ichorCNA.img"

    input:
        tuple val(meta), path(wigfile)
        path gc_wig
        path map_wig
        path pon_file
        path centromere
        path reptime_file

    output:
        tuple val(meta), path("*.cna.seg")           , emit: cna_seg
        tuple val(meta), path("*.seg.txt")           , emit: bins
        tuple val(meta), path("*.seg")      
        tuple val(meta), path("*.correctedDepth.txt"), emit: corrected_depth
        path("*.params.txt")                         , emit: ichorcna_params
        path "**/*genomeWide.pdf"                    , emit: genome_plot
        path "*.RData"
        path "versions.yml"                          , emit: versions

    script:

        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"

        def pon = params.normal_panel ? "--normalPanel ${params.normal_panel}" : ''
        def centro = params.centromere ? "--centromere ${params.centromere}" : ''
        def reptime = params.reptime_file ? "--repTimeWig ${params.reptime_file}": ''
        def ichorcna_script = "/usr/local/bin/runIchorCNA.R"
        def gender = meta.gender ? "--sex ${meta.gender}": ''

        def VERSION = '0.3.2' // WARN: Version information not provided by tool on CLI.

        """
        Rscript ${ichorcna_script} \\
                --id ${prefix} \\
                --WIG $wigfile \\
                --gcWig ${gc_wig} \\
                --mapWig ${map_wig} \\
                $args \\
                ${pon} \\
                ${centro} \\
                ${reptime} \\
                ${gender} \\
                --estimateNormal TRUE \\
                --cores "${task.cpus}" \\
                --outDir ./


        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ichorcna: $VERSION
        END_VERSIONS
        """
}
