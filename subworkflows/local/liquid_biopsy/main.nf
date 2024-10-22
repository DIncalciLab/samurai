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
include { CORRECT_LOGR_ICHORCNA                               } from '../../../modules/local/correct_logR_ichorcna/main'

include { HMMCOPY_READCOUNTER as HMMCOPY_READCOUNTER_ICHORCNA } from '../../../modules/nf-core/hmmcopy/readcounter/main'

workflow LIQUID_BIOPSY {

    take:
        bam_bai // [meta, bam, bai]
        caller

    main:

        ch_versions = Channel.empty()
        ch_reports = Channel.empty()
        ch_bam_bai = bam_bai

        // If we want to build the normal panel
        if (params.build_pon) {

            BUILD_PON(params.pon_path, caller)
            ch_versions = ch_versions.mix(BUILD_PON.out.versions)
            pon_file = BUILD_PON.out.normal_panel

        } else {
            if (!params.normal_panel) {
                if caller == "wisecondorX" ) {
                    error "No PoN specified nor built, but WisecondorX requires it"
                }
                // ichorCNA can work without a PoN, although not optimally
                log.warn "No PoN specified: CNA calling performance can be impacted"
                pon_file = []

            } else {
                pon_file = file(params.normal_panel, checkIfExists: true)
            }
        }

        switch(caller) {
            case "ichorcna":

                gc_wig              = file(params.ichorcna_gc_wig, checkIfExists: true)
                map_wig             = file(params.ichorcna_map_wig, checkIfExists: true)
                centromere          = file(params.ichorcna_centromere_file, checkIfExists: true)
                reptime_file        = params.ichorcna_reptime_wig ? file(params.ichorcna_reptime_wig, checkIfExists: true): []

                // Step 1: Generate coverage wig files
                // To generate Wig Files
                HMMCOPY_READCOUNTER_ICHORCNA(ch_bam_bai)
                // wigfiles = HMMCOPY_READCOUNTER_ICHORCNA.out.wig.map{it -> it[1]}
                ch_versions = ch_versions.mix(HMMCOPY_READCOUNTER_ICHORCNA.out.versions)

                // Step 2: run ichorCNA
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

                // Step 3: produce an aggregate table of the results
                AGGREGATE_ICHORCNA_TABLE (
                    RUN_ICHORCNA.out.ichorcna_params.collect())

                ch_versions = ch_versions.mix(AGGREGATE_ICHORCNA_TABLE.out.versions)
                ch_reports = ch_reports.mix(AGGREGATE_ICHORCNA_TABLE.out.ichorcna_summary)

                RUN_ICHORCNA.out.bins.map{meta, data -> data}
                                     .collectFile(storeDir: "${params.outdir}/ichorcna/",
                                                 name: 'all_segments_ichorcna_gistic.seg',
                                                 keepHeader: true,
                                                 skip: 1)
                                                 .set{gistic_file}

                CORRECT_LOGR_ICHORCNA(gistic_file, AGGREGATE_ICHORCNA_TABLE.out.ichorcna_summary)
                ch_versions = ch_versions.mix(CORRECT_LOGR_ICHORCNA.out.versions)

                corrected_gistic_file = CORRECT_LOGR_ICHORCNA.out.gistic_file
                // Step 4: Aggregate bin-level plots into a single file
                CONCATENATE_BIN_PLOTS(RUN_ICHORCNA.out.genome_plot.collect())
                ch_versions = ch_versions.mix(CONCATENATE_BIN_PLOTS.out.versions)

                break

            case "wisecondorx":
                blacklist = params.wisecondorx_blacklist ? file(params.wisecondorx_blacklist, checkIfExists: true): []

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

                CONVERT_GISTIC_SEG.out.gistic_file
                                    .map{meta, data -> data}
                                    .collectFile(name: "all_segments_wisecondorx_gistic.seg",
                                    skip:1,
                                    storeDir: "${params.outdir}/wisecondorx/" )
                                    .set{ gistic_file }
                ASSEMBLE_WISECONDORX_OUTPUTS(
                            WISECONDORX_PREDICT.out.statistics.collect{
                                meta, result -> result
                            },
                            WISECONDORX_PREDICT.out.calls.collect{
                                meta, result -> result
                            },)
                ch_versions = ch_versions.mix(ASSEMBLE_WISECONDORX_OUTPUTS.out.versions)
                ch_reports = ch_reports.mix(ASSEMBLE_WISECONDORX_OUTPUTS.out.wisecondorx_summary)

                CONVERT_WISECONDORX_IMAGES(
                    WISECONDORX_PREDICT.out.genome_plot.collect{
                        meta, result -> result
                    }
                )
                ch_versions = ch_versions.mix(CONVERT_WISECONDORX_IMAGES.out.versions)

                genome_plot = CONVERT_WISECONDORX_IMAGES.out.genome_plot
                // For compatibility with workflow output
                corrected_gistic_file = gistic_file
                signature_file = gistic_file
                break
            default:
                error "Uknown / unsupported analysis type ${caller}"
        }

    emit:

        normal_panel          = pon_file
        called_segments       = called_segments
        genome_plot           = genome_plot
        summary               = ch_reports
        corrected_gistic_file = corrected_gistic_file
        versions              = ch_versions

}
