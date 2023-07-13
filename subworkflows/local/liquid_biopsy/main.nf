// Import modules
include { CONCATENATE_PDF as CONCATENATE_BIN_PLOTS            } from '../../../modules/local/concatenate_pdf/main'
include { ICHORCNA_PON                                        } from '../../../subworkflows/local/ichorcna_pon/main'
include { RUN_ICHORCNA                                        } from '../../../modules/local/ichorcna/run/main'
include { AGGREGATE_TABLE                                     } from '../../../modules/local/aggregate_table/main'
include { HMMCOPY_READCOUNTER as HMMCOPY_READCOUNTER_ICHORCNA } from '../../../modules/nf-core/hmmcopy/readcounter/main'

workflow LIQUID_BIOPSY {

    take:
        bam_bai
        gc_wig
        map_wig
        centromere
        reptime_file



    main:

        ch_versions = Channel.empty()

        // If we want to build the normal panel
        if (params.build_pon) {
            ICHORCNA_PON(params.pon_path, 
                        gc_wig, map_wig, centromere, reptime_file)
            ch_versions = ch_versions.mix(ICHORCNA_PON.out.versions)
            pon_file = ICHORCNA_PON.out.normal_panel
            } else {
                pon_file = Channel.value(params.normal_panel)
            }

         // To generate Wig Files
        HMMCOPY_READCOUNTER_ICHORCNA(bam_bai)
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

        AGGREGATE_TABLE (
               RUN_ICHORCNA.out.ichorcna_params.collect()
        )

        CONCATENATE_BIN_PLOTS(RUN_ICHORCNA.out.genome_plot.collect())
        ch_versions = ch_versions.mix(CONCATENATE_BIN_PLOTS.out.versions)

    emit:

        normal_panel        = pon_file
        ichorcna_params     = RUN_ICHORCNA.out.ichorcna_params
        called_segments     = RUN_ICHORCNA.out.cna_seg
        genome_plot         = RUN_ICHORCNA.out.genome_plot
        all_bin_plots       = CONCATENATE_BIN_PLOTS.out.genome_plot
        ploidy_summary      = AGGREGATE_TABLE.out.ploidy_summary

        versions            = ch_versions

}
