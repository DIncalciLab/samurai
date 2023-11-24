include { GISTIC2                   } from '../../../modules/local/gistic2/main'
include { ASSEMBLE_GISTIC_OUTPUT        } from '../../../modules/local/assemble_gistic_output/main'

workflow RUN_GISTIC {

    take:
        segmentation_file

    main:
        ch_versions = Channel.empty()
        seg_file = segmentation_file

        GISTIC2(seg_file)
        ch_versions = ch_versions.mix(GISTIC2.out.versions)

        
        ASSEMBLE_GISTIC_OUTPUT(GISTIC2.out.all_lesions)
        ch_versions = ch_versions.mix(ASSEMBLE_GISTIC_OUTPUT.out.versions)
    emit:
        gistic_lesions              = ASSEMBLE_GISTIC_OUTPUT.out.gistic_lesions
        // gistic_genes                = ASSEMBLE_GISTIC_OUTPUT.out.gistic_genes NEED TO FIX THIS
        gistic_log_r                = ASSEMBLE_GISTIC_OUTPUT.out.gistic_log_r
        gistic_cn_states            = ASSEMBLE_GISTIC_OUTPUT.out.gistic_cn_states
        versions                    = ch_versions

}
