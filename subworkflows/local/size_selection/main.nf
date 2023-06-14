// Import local modules
include { BAMPE_FRAGMENTSIZE                              } from '../../../modules/local/deeptools/bampefragmentsize/main'

// Import nf-core modules
include { SAMTOOLS_STATS as SAMTOOLS_STATS_PRE            } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_POST           } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_SIZE_SELECTION } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW                                   } from '../../../modules/nf-core/samtools/view/main'

workflow SIZE_SELECTION {     

    take:
        ch_bam_bai
        fasta
    
    main:
        ch_versions = Channel.empty()

        SAMTOOLS_STATS_PRE(ch_bam_bai, fasta)
        ch_versions = ch_versions.mix(SAMTOOLS_STATS_PRE.out.versions)

        SAMTOOLS_VIEW(ch_bam_bai.map{ it -> [[id:it[1].baseName], it[1], it[2]]},
                            fasta.map{ it -> [[id:it[0].baseName], it] }, []) // size selction 
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)

        SAMTOOLS_INDEX_SIZE_SELECTION(SAMTOOLS_VIEW.out.bam)
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_SIZE_SELECTION.out.versions)

        ch_filtered = SAMTOOLS_VIEW.out.bam.join(SAMTOOLS_INDEX_SIZE_SELECTION.out.bai)

        SAMTOOLS_STATS_POST(ch_filtered, [])
        ch_versions = ch_versions.mix(SAMTOOLS_STATS_POST.out.versions)

        ch_fragment = ch_filtered
            .multiMap {
                it ->
                bam: it[1]
                bai: it[2]
            }

        BAMPE_FRAGMENTSIZE(Channel.value("postfiltering"),
                        ch_fragment.bam.collect(),
                        ch_fragment.bai.collect())
        ch_versions = ch_versions.mix(BAMPE_FRAGMENTSIZE.out.versions)

    emit:
        bam            = SAMTOOLS_VIEW.out.bam
        bai            = SAMTOOLS_INDEX_SIZE_SELECTION.out.bai
        stats_post     = SAMTOOLS_STATS_POST.out.stats
        histogram      = BAMPE_FRAGMENTSIZE.out.histogram 
        size_table     = BAMPE_FRAGMENTSIZE.out.table
        size_raw_table = BAMPE_FRAGMENTSIZE.out.raw_table
        versions       = ch_versions
}