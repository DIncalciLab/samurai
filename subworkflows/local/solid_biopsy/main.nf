// Import modules
include { QDNASEQ                                  } from '../../../modules/local/qdnaseq/main'
include { CONCATENATE_PDF as CONCATENATE_BIN_PLOTS } from '../../../modules/local/concatenate_pdf/main'
include { CONCATENATE_PDF as CONCATENATE_SEG_PLOTS } from '../../../modules/local/concatenate_pdf/main'
include { CONCATENATE_SEG as CONCATENATE_SEGMENTS  } from '../../../modules/local/concatenate_seg/main'
include { CONCATENATE_SEG as CONCATENATE_FILT_SEG  } from '../../../modules/local/concatenate_seg/main'
include { CONCATENATE_SEG as CONCATENATE_SEG_CALLS } from '../../../modules/local/concatenate_seg/main'



// Workfow

workflow SOLID_BIOPSY {

    take:
        ch_bam_bai
        ch_binsize

    main:
        ch_versions = Channel.empty()

        QDNASEQ(ch_bam_bai, 
                ch_binsize)

        ch_versions = ch_versions.mix( QDNASEQ.out.versions.first() )
        
        CONCATENATE_BIN_PLOTS(QDNASEQ.out.bin_plot.collect())
        CONCATENATE_SEG_PLOTS(QDNASEQ.out.segment_plot.collect())
        ch_versions = ch_versions.mix( CONCATENATE_BIN_PLOTS.out.versions.first() )

        CONCATENATE_SEGMENTS(QDNASEQ.out.segments.collect(), "all_segments")
        CONCATENATE_FILT_SEG(QDNASEQ.out.filtered_segments.collect(), "all_filtered_segments")
        CONCATENATE_SEG_CALLS(QDNASEQ.out.called_segments.collect(), "all_called_segments")
        
    emit:
        bins                  = QDNASEQ.out.bins //optional: true
        called_segments       = QDNASEQ.out.called_segments
        segments              = QDNASEQ.out.segments
        processed_data        = QDNASEQ.out.processed_data
        bin_plot              = QDNASEQ.out.bin_plot
        segment_plot          = QDNASEQ.out.segment_plot //optional: true
        
        all_segments          = CONCATENATE_SEGMENTS.out.all_segments
        all_filtered_segments = CONCATENATE_FILT_SEG.out.all_segments
        all_calls             = CONCATENATE_SEG_CALLS.out.all_segments
        all_bin_plots         = CONCATENATE_BIN_PLOTS.out.genome_plot
        all_seg_plots         = CONCATENATE_SEG_PLOTS.out.genome_plot

        versions = ch_versions

}