// Import local modules
include { BAMPE_FRAGMENTSIZE as BAMPE_FRAGMENTSIZE_PRE    } from '../../../modules/local/deeptools/bampefragmentsize/main'
include { BAMPE_FRAGMENTSIZE as BAMPE_FRAGMENTSIZE_POST   } from '../../../modules/local/deeptools/bampefragmentsize/main'

// Import nf-core modules
include { SAMTOOLS_STATS as SAMTOOLS_STATS_PRE            } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_POST           } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_SIZE_SELECTION } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW                                   } from '../../../modules/nf-core/samtools/view/main'

workflow SIZE_SELECTION {

    take:
        ch_bam_bai
        ch_fasta
    main:
        ch_versions = Channel.empty()

        SAMTOOLS_STATS_PRE(ch_bam_bai, ch_fasta)
        ch_versions = ch_versions.mix(SAMTOOLS_STATS_PRE.out.versions)

        bam_bai_pre = ch_bam_bai
            .multiMap {
                it ->
                bam: it[1]
                bai: it[2]
            }

        BAMPE_FRAGMENTSIZE_PRE(Channel.value("prefiltering"),
                        bam_bai_pre.bam.collect(),
                        bam_bai_pre.bai.collect())

        SAMTOOLS_VIEW(ch_bam_bai, ch_fasta, []) // size selection
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)

        SAMTOOLS_INDEX_SIZE_SELECTION(SAMTOOLS_VIEW.out.bam)
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_SIZE_SELECTION.out.versions)

        ch_filtered = SAMTOOLS_VIEW.out.bam.join(SAMTOOLS_INDEX_SIZE_SELECTION.out.bai)

        SAMTOOLS_STATS_POST(ch_filtered, ch_fasta)
        ch_versions = ch_versions.mix(SAMTOOLS_STATS_POST.out.versions)

        ch_fragment = ch_filtered
            .multiMap {
                it ->
                bam: it[1]
                bai: it[2]
            }

        BAMPE_FRAGMENTSIZE_POST(Channel.value("postfiltering"),
                        ch_fragment.bam.collect(),
                        ch_fragment.bai.collect())
        ch_versions = ch_versions.mix(BAMPE_FRAGMENTSIZE_POST.out.versions)

    emit:
        bam            = SAMTOOLS_VIEW.out.bam
        bai            = SAMTOOLS_INDEX_SIZE_SELECTION.out.bai
        stats_pre      = SAMTOOLS_STATS_PRE.out.stats
        stats_post     = SAMTOOLS_STATS_POST.out.stats
        histogram      = BAMPE_FRAGMENTSIZE_POST.out.histogram
        size_table     = BAMPE_FRAGMENTSIZE_POST.out.table
        size_raw_table = BAMPE_FRAGMENTSIZE_POST.out.raw_table
        versions       = ch_versions
}
