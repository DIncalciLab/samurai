include { GISTIC2                   } from '../../../modules/local/gistic2/main'
include { ASSEMBLE_GISTIC_OUTPUT    } from '../../../modules/local/assemble_gistic_output/main'
include { MAFTOOLS                  } from '../../../modules/local/maftools/main'

workflow RUN_GISTIC {

    take:
        gistic_results_dir
        genome // value: genome

    main:
        ch_versions = Channel.empty()

        GISTIC2(gistic_results_dir, genome)
        ch_versions = ch_versions.mix(GISTIC2.out.versions)

        ASSEMBLE_GISTIC_OUTPUT(GISTIC2.out.gistic_results_dir)
        ch_versions = ch_versions.mix(ASSEMBLE_GISTIC_OUTPUT.out.versions)

        MAFTOOLS(GISTIC2.out.all_lesions,
                GISTIC2.out.amplified_genes,
                GISTIC2.out.deleted_genes,
                GISTIC2.out.gistic_score)
        ch_versions = ch_versions.mix(MAFTOOLS.out.versions)

    emit:
        gistic_lesions              = ASSEMBLE_GISTIC_OUTPUT.out.gistic_lesions
        gistic_broad_lesions        = ASSEMBLE_GISTIC_OUTPUT.out.gistic_broad_lesions
        gistic_genes                = ASSEMBLE_GISTIC_OUTPUT.out.gistic_genes
        gistic_log_r                = ASSEMBLE_GISTIC_OUTPUT.out.gistic_log_r
        gistic_cn_states            = ASSEMBLE_GISTIC_OUTPUT.out.gistic_cn_states
        gistic_lesions              = ASSEMBLE_GISTIC_OUTPUT.out.gistic_lesions
        chrom_plot                  = MAFTOOLS.out.chrom_plot
        versions                    = ch_versions

}
