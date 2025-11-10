// Import modules

include { BUILD_PON                    } from '../../../subworkflows/local/build_pon/main'
include { CONVERT_GISTIC_SEG           } from '../../../modules/local/convert_gistic_seg/main'
include { BAM_CNV_WISECONDORX          } from '../../../subworkflows/nf-core/bam_cnv_wisecondorx/main'
include { ICHORCNA                     } from '../../../subworkflows/local/ichorcna/main'
include { ASSEMBLE_WISECONDORX_OUTPUTS } from '../../../modules/local/assemble_wisecondorx_outputs/main'
include { CONVERT_WISECONDORX_IMAGES   } from '../../../modules/local/convert_wisecondorx_images/main'

workflow LIQUID_BIOPSY {
    take:
    ch_bam_bai // [meta, bam, bai]
    caller
    ch_fasta // [meta2, fasta]
    ch_fai // [meta3, fasta]
    ch_normal_panel // channel: normal_panel
    ch_gc_wig // channel: path to GC wig
    ch_map_wig // channel: path to mappability wig
    ch_centromere // channel: path to centromere file
    ch_reptiming // channel: path to reptiming file
    build_pon // bool
    pon_path // value: path
    ch_blacklist // channel: [meta, blacklist]

    main:

    ch_versions = channel.empty()
    ch_reports = channel.empty()

    // If we want to build the normal panel
    if (build_pon) {

        BUILD_PON(pon_path, caller, ch_fasta, ch_fai, ch_gc_wig, ch_map_wig, ch_reptiming, ch_centromere)
        ch_versions = ch_versions.mix(BUILD_PON.out.versions)
        pon_file = BUILD_PON.out.normal_panel.collect()
    }
    else {
        if (!ch_normal_panel) {
            if (caller == "wisecondorx") {
                error("No PoN specified nor built, but WisecondorX requires it")
            }
            else {
                // ichorCNA can work without a PoN, although not optimally
                log.warn("No PoN specified: CNA calling performance can be impacted")
                pon_file = []
            }
        }
        else {
            pon_file = ch_normal_panel.collect()
        }
    }

    if (caller == "ichorcna") {

        ICHORCNA(ch_bam_bai, pon_file, ch_gc_wig, ch_map_wig, ch_centromere, ch_reptiming, ch_fasta)
        ch_versions = ch_versions.mix(ICHORCNA.out.versions)
        ch_reports = ch_reports.mix(ICHORCNA.out.summary)
        called_segments = ICHORCNA.out.ch_segments
        genome_plot = ICHORCNA.out.genome_plot
        corrected_gistic_file = ICHORCNA.out.gistic_file
    }
    else if (caller == "wisecondorx") {
        
        BAM_CNV_WISECONDORX(ch_bam_bai, ch_fasta, ch_fai, pon_file, ch_blacklist)
        ch_versions = ch_versions.mix(BAM_CNV_WISECONDORX.out.versions)

        CONVERT_GISTIC_SEG(
            BAM_CNV_WISECONDORX.out.segments_bed,
            BAM_CNV_WISECONDORX.out.bins_bed,
            BAM_CNV_WISECONDORX.out.aberrations_bed,
        )
        ch_versions = ch_versions.mix(CONVERT_GISTIC_SEG.out.versions)

        called_segments = BAM_CNV_WISECONDORX.out.aberrations_bed
        CONVERT_GISTIC_SEG.out.gistic_file
            .map { _meta, data -> data }
            .collectFile(
                name: "all_segments_wisecondorx_gistic.seg",
                skip: 1,
                storeDir: "${params.outdir}/wisecondorx/",
            )
            .set { gistic_file }

        ASSEMBLE_WISECONDORX_OUTPUTS(
            BAM_CNV_WISECONDORX.out.chr_statistics.collect { _meta, result ->
                result
            },
            BAM_CNV_WISECONDORX.out.aberrations_bed.collect { _meta, result ->
                result
            },
        )
        ch_versions = ch_versions.mix(ASSEMBLE_WISECONDORX_OUTPUTS.out.versions)
        ch_reports = ch_reports.mix(ASSEMBLE_WISECONDORX_OUTPUTS.out.wisecondorx_summary)

        CONVERT_WISECONDORX_IMAGES(
            BAM_CNV_WISECONDORX.out.genome_plot.collect { _meta, result ->
                result
            }
        )
        ch_versions = ch_versions.mix(CONVERT_WISECONDORX_IMAGES.out.versions)

        genome_plot = CONVERT_WISECONDORX_IMAGES.out.genome_plot
        // For compatibility with workflow output
        corrected_gistic_file = gistic_file
    }
    else {
        error("Uknown / unsupported analysis type ${caller}")
    }

    emit:
    normal_panel          = pon_file
    called_segments       = called_segments
    genome_plot           = genome_plot
    summary               = ch_reports
    corrected_gistic_file = corrected_gistic_file
    versions              = ch_versions
}
