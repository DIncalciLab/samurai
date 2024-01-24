// Import modules
include { QDNASEQ                                               } from '../../../modules/local/qdnaseq/main'
include { ASCAT_SC                                              } from '../../../modules/local/ascat_sc/main'
include { CONCATENATE_PDF as CONCATENATE_QDNASEQ_PLOTS          } from '../../../modules/local/concatenate_pdf/main'
include { CONCATENATE_PDF as CONCATENATE_ASCATSC_PLOTS          } from '../../../modules/local/concatenate_pdf/main'
include { CONCATENATE_PDF as CONCATENATE_ASCATSC_REFITTED_PLOTS } from '../../../modules/local/concatenate_pdf/main'
include { CREATE_QDNASEQ_SUMMARY                                } from '../../../modules/local/create_qdnaseq_summary/main'
include { CREATE_ASCATSC_SUMMARY                                } from '../../../modules/local/create_ascatsc_summary/main'
include { CIN_SIGNATURE_QUANTIFICATION                          } from '../../../modules/local/cin_signature_quantification/main'

// Workfow

workflow SOLID_BIOPSY {

    take:
        ch_bam_bai // channel
        caller // value: either "qdnaseq" or "ascat_sc"
        binsize // value, bin size in bp
        genome // value, genome to use

    main:
        ch_versions = Channel.empty()
        ch_reports = Channel.empty()

        switch(caller) {
            case "qdnaseq":
                QDNASEQ(ch_bam_bai, binsize, genome)
                ch_versions = ch_versions.mix(QDNASEQ.out.versions.first())
                CONCATENATE_QDNASEQ_PLOTS(QDNASEQ.out.segment_plot.collect())
                ch_versions = ch_versions.mix(CONCATENATE_QDNASEQ_PLOTS.out.versions)
                all_seg_plot = CONCATENATE_QDNASEQ_PLOTS.out.genome_plot

                QDNASEQ.out.segments
                            .collectFile(storeDir: "${params.outdir}/qdnaseq/",
                                         name: 'all_segments.seg',
                                         keepHeader: true,
                                         skip: 1)
                            .set{ ch_segments }
                QDNASEQ.out.summary_table
                            .collectFile(storeDir: "${params.outdir}/qdnaseq/",
                                         name: 'qdnaseq_summary_mqc.txt',
                                         keepHeader: true,
                                         skip: 1)
                            .set{ qdnaseq_summary }
                //CREATE_QDNASEQ_SUMMARY(qdnaseq_summary)
                //ch_versions = ch_versions.mix(CREATE_QDNASEQ_SUMMARY.out.versions)
                ch_reports = ch_reports.mix(qdnaseq_summary)
                //TODO: Generate the GISTIC output here? If not, remove the step from the module
                break
            case "ascat_sc":
                ASCAT_SC(ch_bam_bai, binsize, genome)
                ch_versions = ch_versions.mix(ASCAT_SC.out.versions)

                if (params.ascat_sc_predict_refit) {
                    CONCATENATE_ASCATSC_REFITTED_PLOTS(ASCAT_SC.out.profiles_refitted.collect())
                    ch_versions = ch_versions.mix(CONCATENATE_ASCATSC_REFITTED_PLOTS.out.versions)
                    all_seg_plot = CONCATENATE_ASCATSC_REFITTED_PLOTS.out.genome_plot

                } else {
                    CONCATENATE_ASCATSC_PLOTS(ASCAT_SC.out.profiles_plot.collect())
                    ch_versions = ch_versions.mix(CONCATENATE_ASCATSC_PLOTS.out.versions)
                    all_seg_plot = CONCATENATE_ASCATSC_PLOTS.out.genome_plot}

                ASCAT_SC.out.segments
                            .collectFile(storeDir: "${params.outdir}/ascat_sc/",
                                                      name: 'all_segments.seg',
                                                      keepHeader: true,
                                                      skip: 1)
                            .set{ ch_segments }

                ASCAT_SC.out.summary_table
                            .collectFile(storeDir: "${params.outdir}/ascat_sc/",
                                 name: 'ascatsc_summary_mqc.txt',
                                 keepHeader: true,
                                 skip: 1)
                            .set { ascatsc_summary }

                ASCAT_SC.out.sig_file
                            .collectFile(storeDir: "${params.outdir}/ascat_sc/",
                                    name: 'all_segments_ascat_sc_signatures.seg',
                                    keepHeader: true,
                                    skip: 1)
                            .set{signature_file}

                if (params.compute_signatures) {
                        CIN_SIGNATURE_QUANTIFICATION(signature_file)
                        ch_versions = ch_versions.mix(CIN_SIGNATURE_QUANTIFICATION.out.versions)
                }

                //CREATE_ASCATSC_SUMMARY(ascatsc_summary)
                ch_reports = ch_reports.mix(ascatsc_summary)

                //ASCAT_SC.out.gistic_file
                //            .collectFile(storeDir: "${params.outdir}/ascat_sc/",
                //                    name: 'all_segments_ascat_sc_gistic.seg',
                //                    keepHeader: true,
                //                    skip: 1)
                //            .set{gistic_file}
                break
            default:
                error "Unknown CNV caller ${caller}"


        }

    emit:
        ch_segments     = ch_segments
        summary         = ch_reports
        versions        = ch_versions
}
