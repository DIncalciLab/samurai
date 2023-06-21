// Import modules
include { QDNASEQ                                  } from '../../../modules/local/qdnaseq/main'
include { CONCATENATE_PDF as CONCATENATE_BIN_PLOTS } from '../../../modules/local/concatenate_pdf/main'
include { CONCATENATE_PDF as CONCATENATE_SEG_PLOTS } from '../../../modules/local/concatenate_pdf/main'




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

        QDNASEQ.out.segments
                    .collectFile(storeDir: "${params.outdir}/qdnaseq/", 
                                 name: 'all_segments.seg', 
                                 keepHeader: true, 
                                 skip: 1)
                                .set{ all_seg_ch}


        QDNASEQ.out.called_segments
                    .collectFile(storeDir: "${params.outdir}/qdnaseq/",
                                 name: 'all_called_segments.seg', 
                                 keepHeader: true, 
                                 skip: 1)
                                 .set{ all_called_seg_ch}
        QDNASEQ.out.summary_table
                    .collectFile(storeDir: "${params.outdir}/qdnaseq/",
                                 name: 'summary_qdnaseq.txt', 
                                 keepHeader: true, 
                                 skip: 1)
                                 .set{ summary_all_samples}

        
    emit:
        bins                  = QDNASEQ.out.bins //optional: true
        called_segments       = QDNASEQ.out.called_segments
        segments              = QDNASEQ.out.segments
        processed_data        = QDNASEQ.out.processed_data
        bin_plot              = QDNASEQ.out.bin_plot
        segment_plot          = QDNASEQ.out.segment_plot //optional: true
        
        all_segments          = all_seg_ch
        all_calls             = all_called_seg_ch
        summary               = summary_all_samples
        all_bin_plots         = CONCATENATE_BIN_PLOTS.out.genome_plot
        all_seg_plots         = CONCATENATE_SEG_PLOTS.out.genome_plot

        versions = ch_versions

}