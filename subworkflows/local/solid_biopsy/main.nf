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
        ch_bam_bai

    main:
        ch_versions = Channel.empty()

        //TODO: Use a single parameter, which goes for either QDNAseq or ASCAT.sc
        // Something like '--analysis_method {qdnaseq,ascatsc}'
        if (params.qdnaseq) {

            binfile             = file(params.binfile)

            QDNASEQ(ch_bam_bai)
            ch_versions = ch_versions.mix(QDNASEQ.out.versions.first() )

            CONCATENATE_QDNASEQ_PLOTS(QDNASEQ.out.bin_plot.collect())
            ch_versions = ch_versions.mix(CONCATENATE_QDNASEQ_PLOTS.out.versions)

            all_seg_plot = CONCATENATE_QDNASEQ_PLOTS.out.genome_plot

            QDNASEQ.out.segments
                        .collectFile(storeDir: "${params.outdir}/qdnaseq/",
                                    name: 'all_segments.seg',
                                    keepHeader: true,
                                    skip: 1)
                                    .set{ all_seg_ch }

            QDNASEQ.out.summary_table
                        .collectFile(storeDir: "${params.outdir}/qdnaseq/",
                                    name: 'qdnaseq_summary.txt',
                                    keepHeader: true,
                                    skip: 1)
                                    .set{ qdnaseq_summary }

            CREATE_QDNASEQ_SUMMARY(qdnaseq_summary)
            ch_versions = ch_versions.mix(CREATE_QDNASEQ_SUMMARY.out.versions)

            summary_multiqc = CREATE_QDNASEQ_SUMMARY.out.qdnaseq_summary
        }

        if (params.ascat_sc) {

            ASCAT_SC(ch_bam_bai)
            ch_versions = ch_versions.mix(ASCAT_SC.out.versions )

            CONCATENATE_ASCATSC_PLOTS(ASCAT_SC.out.profiles_plot.collect())
            ch_versions = ch_versions.mix(CONCATENATE_ASCATSC_PLOTS.out.versions)
            all_seg_plot = CONCATENATE_ASCATSC_PLOTS.out.genome_plot

            ASCAT_SC.out.segments
                        .collectFile(storeDir: "${params.outdir}/ascat_sc/",
                                    name: 'all_segments.seg',
                                    keepHeader: true,
                                    skip: 1)
                                    .set{ all_seg_ch }

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
            ch_versions = ch_versions.mix(CREATE_ASCATSC_SUMMARY.out.versions)

            summary_multiqc = CREATE_ASCATSC_SUMMARY.out.ascatsc_summary

            if (params.quantify_signatures) {
                QUANTIFY_CIN_SIGNATURES(signature_file)
                ch_versions = ch_versions.mix( QUANTIFY_CIN_SIGNATURES.out.versions)
            }



        }

    //TODO: Differentiate between ASCAT.sc and QDNAseq outputs somehow
    // Otherwise, what happens if you run *both* ?
    emit:
        all_seg_ch      = all_seg_ch
        summary         = summary_multiqc
        versions        = ch_versions
}
