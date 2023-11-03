// Import modules
include { QDNASEQ                                      } from '../../../modules/local/qdnaseq/main'
include { ASCAT_SC                                     } from '../../../modules/local/ascat_sc/main'
include { CONCATENATE_PDF as CONCATENATE_QDNASEQ_PLOTS } from '../../../modules/local/concatenate_pdf/main'
include { CONCATENATE_PDF as CONCATENATE_ASCATSC_PLOTS } from '../../../modules/local/concatenate_pdf/main'
include { CREATE_QDNASEQ_SUMMARY                       } from '../../../modules/local/create_qdnaseq_summary/main'
include { CREATE_ASCATSC_SUMMARY                       } from '../../../modules/local/create_ascatsc_summary/main'
include { QUANTIFY_CIN_SIGNATURES                      } from '../../../modules/local/quantify_cin_signatures/main'

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
                                         name: 'qdnaseq_summary.txt',
                                         keepHeader: true,
                                         skip: 1)
                            .set{ qdnaseq_summary }
                CREATE_QDNASEQ_SUMMARY(qdnaseq_summary)
                ch_versions = ch_versions.mix(CREATE_QDNASEQ_SUMMARY.out.versions)
                ch_reports = ch_reports.mix(CREATE_QDNASEQ_SUMMARY.out.summary)
                //TODO: Generate the GISTIC output here? If not, remove the step from the module
                break
            case "ascat_sc":
                ASCAT_SC(ch_bam_bai, binsize, genome)
                ch_versions = ch_versions.mix(ASCAT_SC.out.versions)

                CONCATENATE_ASCATSC_PLOTS(ASCAT_SC.out.profiles_plot.collect())
                ch_versions = ch_versions.mix(CONCATENATE_ASCATSC_PLOTS.out.versions)
                all_seg_plot = CONCATENATE_ASCATSC_PLOTS.out.genome_plot

                ASCAT_SC.out.segments
                    .collectFile(storeDir: "${params.outdir}/ascat_sc/",
                                              name: 'all_segments.seg',
                                              keepHeader: true,
                                              skip: 1)
                    .set{ ch_segments }

                ASCAT_SC.out.summary_table
                    .collectFile(storeDir: "${params.outdir}/ascat_sc/",
                                 name: 'ascatsc_summary.txt',
                                 keepHeader: true,
                                 skip: 1)
                    .set { ascatsc_summary }

                ASCAT_SC.out.sig_file
                        .collectFile(storeDir: "${params.outdir}/ascat_sc/",
                                    name: 'segments_sig_extraction.seg',
                                    keepHeader: true,
                                    skip: 1)
                                    .set{signature_file}

                CREATE_ASCATSC_SUMMARY(ascatsc_summary)
                ch_reports = ch_reports.mix(CREATE_ASCATSC_SUMMARY.out.summary)
                //TODO: Generate the GISTIC output here? If not, remove the step from the module
                break
            default:
                error "Unknown CNV caller ${caller}"
            if (params.quantify_signatures) {
                QUANTIFY_CIN_SIGNATURES(signature_file)
                ch_versions = ch_versions.mix( QUANTIFY_CIN_SIGNATURES.out.versions)
            }


        }

    emit:
        ch_segments     = ch_segments
        summary         = ch_reports
        versions        = ch_versions
}
