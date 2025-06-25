include { RUN_ICHORCNA                                        } from '../../../modules/local/ichorcna/run/main'
include { AGGREGATE_ICHORCNA_TABLE                            } from '../../../modules/local/aggregate_ichorcna_table/main'
include { HMMCOPY_READCOUNTER as HMMCOPY_READCOUNTER_ICHORCNA } from '../../../modules/nf-core/hmmcopy/readcounter/main'
include { CORRECT_LOGR_ICHORCNA                               } from '../../../modules/local/correct_logR_ichorcna/main'
include { CONCATENATE_PDF as CONCATENATE_BIN_PLOTS            } from '../../../modules/local/concatenate_pdf/main'

workflow ICHORCNA {
    take:
    ch_bam_bai      // [meta, bam, bai]
    ch_normal_panel // Channel: path (optional)
    ch_gc_wig       // Channel: path
    ch_map_wig      // Channel: path
    ch_centromere   // Channel: path
    ch_reptime_wig  // Channel: path (optional)

    main:

    ch_versions = Channel.empty()
    ch_reports = Channel.empty()
    // Step 1: Generate coverage wig files
    // To generate Wig Files
    HMMCOPY_READCOUNTER_ICHORCNA(ch_bam_bai)
    ch_versions = ch_versions.mix(HMMCOPY_READCOUNTER_ICHORCNA.out.versions)
    // Step 2: run ichorCNA
    RUN_ICHORCNA(
        HMMCOPY_READCOUNTER_ICHORCNA.out.wig,
        ch_gc_wig,
        ch_map_wig,
        ch_normal_panel,
        ch_centromere,
        ch_reptime_wig,
    )
    ch_versions = ch_versions.mix(RUN_ICHORCNA.out.versions)

    called_segments = RUN_ICHORCNA.out.cna_seg
    bins = RUN_ICHORCNA.out.bins
    genome_plot = RUN_ICHORCNA.out.genome_plot

    // Step 3: produce an aggregate table of the results
    AGGREGATE_ICHORCNA_TABLE(
        RUN_ICHORCNA.out.ichorcna_params.collect()
    )

    ch_versions = ch_versions.mix(AGGREGATE_ICHORCNA_TABLE.out.versions)
    ch_reports = ch_reports.mix(AGGREGATE_ICHORCNA_TABLE.out.ichorcna_summary)

    RUN_ICHORCNA.out.cna_seg
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
    // Step 4: Aggregate bin-level plots into a single file
    CONCATENATE_BIN_PLOTS(RUN_ICHORCNA.out.genome_plot.collect())
    ch_versions = ch_versions.mix(CONCATENATE_BIN_PLOTS.out.versions)

    emit:
    versions    = ch_versions
    summary     = ch_reports
    ch_segments = called_segments
    ch_bins     = bins
    gistic_file = corrected_gistic_file
    genome_plot = genome_plot
}
