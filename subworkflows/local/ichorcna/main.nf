include { ICHORCNA_RUN                                        } from '../../../modules/nf-core/ichorcna/run/main'
include { AGGREGATE_ICHORCNA_TABLE                            } from '../../../modules/local/aggregate_ichorcna_table/main'
include { HMMCOPY_READCOUNTER as HMMCOPY_READCOUNTER_ICHORCNA } from '../../../modules/nf-core/hmmcopy/readcounter/main'
include { CORRECT_LOGR_ICHORCNA                               } from '../../../modules/local/correct_logR_ichorcna/main'
include { CONCATENATE_PDF as CONCATENATE_BIN_PLOTS            } from '../../../modules/local/concatenate_pdf/main'
include { PLOT_ICHORCNA                                       } from '../../../modules/local/plot_ichorcna/main'

workflow ICHORCNA {
    take:
    ch_bam_bai // [meta, bam, bai]
    ch_normal_panel // Channel: path (optional)
    ch_gc_wig // Channel: path
    ch_map_wig // Channel: path
    ch_centromere // Channel: path
    ch_reptime_wig // Channel: path (optional)
    ch_fasta // Channel: [meta, fasta]

    main:

    ch_versions = channel.empty()
    ch_reports = channel.empty()
    // Step 1: Generate coverage wig files
    // To generate Wig Files
    HMMCOPY_READCOUNTER_ICHORCNA(ch_bam_bai, ch_fasta)
    ch_versions = ch_versions.mix(HMMCOPY_READCOUNTER_ICHORCNA.out.versions)
    // Step 2: run ichorCNA

    ICHORCNA_RUN(
        HMMCOPY_READCOUNTER_ICHORCNA.out.wig,
        ch_gc_wig,
        ch_map_wig,
        [],
        ch_normal_panel,
        ch_centromere,
        ch_reptime_wig,
        [],
    )

    ch_versions = ch_versions.mix(ICHORCNA_RUN.out.versions)

    called_segments = ICHORCNA_RUN.out.seg_txt
    bins = ICHORCNA_RUN.out.cna_seg
    genome_plot = ICHORCNA_RUN.out.genome_plot

    // Step 3: produce an aggregate table of the results
    AGGREGATE_ICHORCNA_TABLE(
        ICHORCNA_RUN.out.ichorcna_params.collect { _meta, params -> params }
    )

    ch_versions = ch_versions.mix(AGGREGATE_ICHORCNA_TABLE.out.versions)
    ch_reports = ch_reports.mix(AGGREGATE_ICHORCNA_TABLE.out.ichorcna_summary)

    ICHORCNA_RUN.out.cna_seg
        .map { _meta, data -> data }
        .collectFile(
            storeDir: "${params.outdir}/ichorcna/",
            name: 'all_segments_ichorcna_gistic.seg',
            keepHeader: true,
            skip: 1,
        )
        .set { gistic_file }

    CORRECT_LOGR_ICHORCNA(gistic_file, AGGREGATE_ICHORCNA_TABLE.out.ichorcna_summary)
    ch_versions = ch_versions.mix(CORRECT_LOGR_ICHORCNA.out.versions)

    corrected_gistic_file = CORRECT_LOGR_ICHORCNA.out.gistic_file
    PLOT_ICHORCNA(
        ICHORCNA_RUN.out.seg_txt,
        ICHORCNA_RUN.out.cna_seg,
        ICHORCNA_RUN.out.ichorcna_params.map { _meta, param -> param },
    )
    ch_versions = ch_versions.mix(PLOT_ICHORCNA.out.versions)

    // Step 4: Aggregate bin-level plots into a single file
    CONCATENATE_BIN_PLOTS(ICHORCNA_RUN.out.genome_plot.collect { _meta, plot -> plot })
    ch_versions = ch_versions.mix(CONCATENATE_BIN_PLOTS.out.versions)

    emit:
    versions    = ch_versions
    summary     = ch_reports
    ch_segments = called_segments
    ch_bins     = bins
    gistic_file = corrected_gistic_file
    genome_plot = genome_plot
}
