process RUN_ICHORCNA {

    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    container "quay.io/dincalcilab/ichorcna:0.4.0-2ab0be2"

    input:
        tuple val(meta), path(wigfile)
        path gc_wig
        path map_wig
        path pon_file
        path centromere
        path reptime_file

    output:
        tuple val(meta), path("*.cna.seg")           , emit: bins
        tuple val(meta), path("*.seg.txt")           , emit: cna_seg
        tuple val(meta), path("*.seg")
        tuple val(meta), path("*.correctedDepth.txt"), emit: corrected_depth
        path("*.params.txt")                         , emit: ichorcna_params
        path "**/*genomeWide.pdf"                    , emit: genome_plot
        path "**/*_genomeWide_all_sols.pdf"
        path "*.RData"
        path "versions.yml"                          , emit: versions

    script:

        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"

        def pon = pon_file ? "--normalPanel ${pon_file}" : ''
        def centro = params.ichorcna_centromere_file ? "--centromere ${params.ichorcna_centromere_file}" : ''
        def reptime = params.ichorcna_reptime_wig ? "--repTimeWig ${params.ichorcna_reptime_wig}": ''
        def ichorcna_script = "runIchorCNA.R"
        def gender = meta.gender ? "--sex ${meta.gender}": ''

        def VERSION = '0.3.2' // WARN: Version information not provided by tool on CLI.

        """
        ${ichorcna_script} \\
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
