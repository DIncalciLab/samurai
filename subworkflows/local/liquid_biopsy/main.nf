// Import modules

include { BUILD_PON                                           } from '../../../subworkflows/local/build_pon/main'
include { RUN_ICHORCNA                                        } from '../../../modules/local/ichorcna/run/main'
include { AGGREGATE_ICHORCNA_TABLE                            } from '../../../modules/local/aggregate_ichorcna_table/main'
include { WISECONDORX_CONVERT                                 } from '../../../modules/local/wisecondorx/convert/main'
include { WISECONDORX_PREDICT                                 } from '../../../modules/local/wisecondorx/predict/main'
include { CONVERT_GISTIC_SEG                                  } from '../../../modules/local/convert_gistic_seg/main'
include { ASSEMBLE_WISECONDORX_OUTPUTS                        } from '../../../modules/local/assemble_wisecondorx_outputs/main'
include { CONVERT_WISECONDORX_IMAGES                          } from '../../../modules/local/convert_wisecondorx_images/main'
include { CONCATENATE_PDF as CONCATENATE_BIN_PLOTS            } from '../../../modules/local/concatenate_pdf/main'

include { HMMCOPY_READCOUNTER as HMMCOPY_READCOUNTER_ICHORCNA } from '../../../modules/nf-core/hmmcopy/readcounter/main'

workflow LIQUID_BIOPSY {

    take:
        bam_bai
        analysis_type

    main:

        ch_versions = Channel.empty()
        ch_bam_bai = bam_bai

        // If we want to build the normal panel
        if (params.build_pon) {

            BUILD_PON(params.pon_path)
            ch_versions = ch_versions.mix(BUILD_PON.out.versions)
            pon_file = BUILD_PON.out.normal_panel

            } else {
                pon_file = Channel.value(params.normal_panel)
            }

        switch(analysis_type) {
            case "liquid_biopsy_ichorcna":

                gc_wig              = Channel.value(params.gc_wig)
                map_wig             = Channel.value(params.map_wig)
                centromere          = Channel.value(params.centromere)
                reptime_file        = Channel.value(params.reptime_file)

                // Step 1: Generate coverage wig files
                // To generate Wig Files
                HMMCOPY_READCOUNTER_ICHORCNA(ch_bam_bai)
                wigfiles = HMMCOPY_READCOUNTER_ICHORCNA.out.wig.map{it -> it[1]}
                ch_versions = ch_versions.mix(HMMCOPY_READCOUNTER_ICHORCNA.out.versions)
                RUN_ICHORCNA (
                    HMMCOPY_READCOUNTER_ICHORCNA.out.wig,
                    gc_wig,
                    map_wig,
                    pon_file,
                    centromere,
                    reptime_file
                )
                ch_versions = ch_versions.mix(RUN_ICHORCNA.out.versions)

                called_segments     = RUN_ICHORCNA.out.cna_seg
                genome_plot         = RUN_ICHORCNA.out.genome_plot

                AGGREGATE_ICHORCNA_TABLE (
                    RUN_ICHORCNA.out.ichorcna_params.collect())
                ch_versions = ch_versions.mix(AGGREGATE_ICHORCNA_TABLE.out.versions)

                summary = AGGREGATE_ICHORCNA_TABLE.out.ichorcna_summary

                CONCATENATE_BIN_PLOTS(RUN_ICHORCNA.out.genome_plot.collect())
                ch_versions = ch_versions.mix(CONCATENATE_BIN_PLOTS.out.versions)
                break

            case "wisecondorx":
                blacklist = params.blacklist ? file(params.blacklist, checkIfExists: true): []

                WISECONDORX_CONVERT(ch_bam_bai)
                ch_versions = ch_versions.mix(WISECONDORX_CONVERT.out.versions)
                ch_npz = WISECONDORX_CONVERT.out.npz

                WISECONDORX_PREDICT(ch_npz, blacklist, pon_file)
                ch_versions = ch_versions.mix(WISECONDORX_PREDICT.out.versions)

                CONVERT_GISTIC_SEG(WISECONDORX_PREDICT.out.segments,
                        WISECONDORX_PREDICT.out.bins,
                        WISECONDORX_PREDICT.out.calls)
                ch_versions = ch_versions.mix(CONVERT_GISTIC_SEG.out.versions)

                called_segments     = WISECONDORX_PREDICT.out.calls

                CONVERT_GISTIC_SEG.out.segfile
                                    .map{meta, data -> data}
                                    .collectFile(name: "all_segments.seg",
                                                    skip: 1,
                                                    storeDir: "${params.outdir}" )
                ASSEMBLE_WISECONDORX_OUTPUTS(
                            WISECONDORX_PREDICT.out.statistics.collect{meta, result -> result},
                            WISECONDORX_PREDICT.out.calls.collect{meta, result -> result},)
                ch_versions = ch_versions.mix(ASSEMBLE_WISECONDORX_OUTPUTS.out.versions)

                summary = ASSEMBLE_WISECONDORX_OUTPUTS.out.wisecondorx_summary

                CONVERT_WISECONDORX_IMAGES(WISECONDORX_PREDICT.out.genome_plot.collect{meta, result -> result})
                ch_versions = ch_versions.mix(CONVERT_WISECONDORX_IMAGES.out.versions)

                genome_plot         = CONVERT_WISECONDORX_IMAGES.out.genome_plot
                break
            default:
                error "Uknown / unsupported analysis type ${analysis_type}"
        }

    emit:

        normal_panel        = pon_file
        called_segments     = called_segments
        genome_plot         = genome_plot
        summary             = summary

        versions            = ch_versions

}
