// Import modules
include { QDNASEQ                                               } from '../../../modules/local/qdnaseq/main'
include { ASCAT_SC                                              } from '../../../modules/local/ascat_sc/main'
include { CONCATENATE_PDF as CONCATENATE_QDNASEQ_PLOTS          } from '../../../modules/local/concatenate_pdf/main'
include { CONCATENATE_PDF as CONCATENATE_ASCATSC_PLOTS          } from '../../../modules/local/concatenate_pdf/main'
include { CONCATENATE_PDF as CONCATENATE_ASCATSC_REFITTED_PLOTS } from '../../../modules/local/concatenate_pdf/main'
include { CREATE_QDNASEQ_SUMMARY                                } from '../../../modules/local/create_qdnaseq_summary/main'
include { CREATE_ASCATSC_SUMMARY                                } from '../../../modules/local/create_ascatsc_summary/main'
include { CIN_SIGNATURE_QUANTIFICATION                          } from '../../../modules/local/cin_signature_quantification/main'
include { HRDCNA                                                } from '../../../modules/local/hrdcna/main'
include { BUILD_PON                                             } from '../../../subworkflows/local/build_pon/main'
include { RUN_ICHORCNA                                          } from '../../../modules/local/ichorcna/run/main'
include { AGGREGATE_ICHORCNA_TABLE                              } from '../../../modules/local/aggregate_ichorcna_table/main'
include { HMMCOPY_READCOUNTER as HMMCOPY_READCOUNTER_ICHORCNA   } from '../../../modules/nf-core/hmmcopy/readcounter/main'
include { CORRECT_LOGR_ICHORCNA                                 } from '../../../modules/local/correct_logR_ichorcna/main'
include { CONCATENATE_PDF as CONCATENATE_BIN_PLOTS              } from '../../../modules/local/concatenate_pdf/main'
include { ICHORCNA                                              } from "../ichorcna/main.nf"
// Workfow

workflow SOLID_BIOPSY {
    take:
    ch_bam_bai // channel
    caller     // value: either "qdnaseq" or "ascat_sc"
    binsize    // value, bin size in bp
    genome     // value, genome to use

    main:
    ch_versions = Channel.empty()
    ch_reports = Channel.empty()

    if (caller == "qdnaseq") {
        QDNASEQ(ch_bam_bai, binsize, genome)
        ch_versions = ch_versions.mix(QDNASEQ.out.versions.first())
        CONCATENATE_QDNASEQ_PLOTS(QDNASEQ.out.segment_plot.collect())
        ch_versions = ch_versions.mix(CONCATENATE_QDNASEQ_PLOTS.out.versions)
        all_seg_plot = CONCATENATE_QDNASEQ_PLOTS.out.genome_plot

        QDNASEQ.out.segments
            .collectFile(
                storeDir: "${params.outdir}/qdnaseq/",
                name: 'all_segments.seg',
                keepHeader: true,
                skip: 1,
            )
            .set { ch_segments }
        QDNASEQ.out.summary_table
            .collectFile(
                storeDir: "${params.outdir}/qdnaseq/",
                name: 'qdnaseq_summary_mqc.txt',
                keepHeader: true,
                skip: 1,
            )
            .set { qdnaseq_summary }
        //CREATE_QDNASEQ_SUMMARY(qdnaseq_summary)
        //ch_versions = ch_versions.mix(CREATE_QDNASEQ_SUMMARY.out.versions)
        ch_reports = ch_reports.mix(qdnaseq_summary)
        corrected_gistic_file = ch_segments
    }
    else if (caller == "ascat_sc") {
        ASCAT_SC(ch_bam_bai, binsize, genome)
        ch_versions = ch_versions.mix(ASCAT_SC.out.versions)

        if (params.ascat_sc_predict_refit) {
            CONCATENATE_ASCATSC_REFITTED_PLOTS(ASCAT_SC.out.profiles_refitted.collect())
            ch_versions = ch_versions.mix(CONCATENATE_ASCATSC_REFITTED_PLOTS.out.versions)
            all_seg_plot = CONCATENATE_ASCATSC_REFITTED_PLOTS.out.genome_plot
        }
        else {
            CONCATENATE_ASCATSC_PLOTS(ASCAT_SC.out.profiles_plot.collect())
            ch_versions = ch_versions.mix(CONCATENATE_ASCATSC_PLOTS.out.versions)
            all_seg_plot = CONCATENATE_ASCATSC_PLOTS.out.genome_plot
        }

        ASCAT_SC.out.segments
            .collectFile(
                storeDir: "${params.outdir}/ascat_sc/",
                name: 'all_segments.seg',
                keepHeader: true,
                skip: 1,
            )
            .set { ch_segments }

        ASCAT_SC.out.summary_table
            .collectFile(
                storeDir: "${params.outdir}/ascat_sc/",
                name: 'ascatsc_summary_mqc.txt',
                keepHeader: true,
                skip: 1,
            )
            .set { ascatsc_summary }

        ASCAT_SC.out.sig_file
            .collectFile(
                storeDir: "${params.outdir}/ascat_sc/",
                name: 'all_segments_ascat_sc_signatures.seg',
                keepHeader: true,
                skip: 1,
            )
            .set { signature_file }

        if (params.compute_signatures) {
            CIN_SIGNATURE_QUANTIFICATION(signature_file)
            ch_versions = ch_versions.mix(CIN_SIGNATURE_QUANTIFICATION.out.versions)
            ch_reports = ch_reports.mix(CIN_SIGNATURE_QUANTIFICATION.out.sig_activity_plot)
        }

        if (params.hrdcna_compute_score) {
            HRDCNA(signature_file)
            ch_versions = ch_versions.mix(HRDCNA.out.versions)
            ch_reports = ch_reports.mix(HRDCNA.out.hrdcna_summary)
        }

        //CREATE_ASCATSC_SUMMARY(ascatsc_summary)
        ch_reports = ch_reports.mix(ascatsc_summary)

        ASCAT_SC.out.gistic_file
            .collectFile(
                storeDir: "${params.outdir}/ascat_sc/",
                name: 'all_segments_ascat_sc_gistic.seg',
                keepHeader: true,
                skip: 1,
            )
            .set { corrected_gistic_file }
    }
    else if (caller == "ichorcna") {

        //FIXME: Massive duplication here with liquid biopsy,
        // we may want to see if we can use sub-sub workflows

        // If we want to build the normal panel
        if (params.build_pon) {

            BUILD_PON(params.pon_path, caller)
            ch_versions = ch_versions.mix(BUILD_PON.out.versions)
            pon_file = BUILD_PON.out.normal_panel
        }
        else {
            if (!params.normal_panel) {
                // ichorCNA can work without a PoN, although not optimally
                log.warn("No PoN specified: CNA calling performance may be impacted")
                pon_file = []
            }
            else {
                pon_file = file(params.normal_panel, checkIfExists: true)
            }
        }

        // FIXME: We shouldn't rely on params here
        gc_wig = file(params.ichorcna_gc_wig, checkIfExists: true)
        map_wig = file(params.ichorcna_map_wig, checkIfExists: true)
        centromere = file(params.ichorcna_centromere_file, checkIfExists: true)
        reptime_file = params.ichorcna_reptime_wig ? file(params.ichorcna_reptime_wig, checkIfExists: true) : []

        ICHORCNA(
            ch_bam_bai,
            pon_file,
            gc_wig,
            map_wig,
            centromere,
            reptime_file,
        )

        ch_segments = ICHORCNA.out.ch_segments
        all_seg_plot = ICHORCNA.out.genome_plot
        corrected_gistic_file = ICHORCNA.out.gistic_file
        ch_reports = ch_versions.mix(ICHORCNA.out.summary)
        ch_versions = ch_versions.mix(ICHORCNA.out.versions)

    }
    else {
        error("Unknown CNV caller ${caller}")
    }

    emit:
    ch_segments = ch_segments
    summary     = ch_reports
    gistic_file = corrected_gistic_file
    versions    = ch_versions
    all_seg_plot = all_seg_plot
}
