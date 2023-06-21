// Import nf-core modules

include { BWA_INDEX                   } from '../../../modules/nf-core/bwa/index/main'
include { BWAMEM2_INDEX               } from '../../../modules/nf-core/bwamem2/index/main'
include { SAMTOOLS_FAIDX              } from '../../../modules/nf-core/samtools/faidx/main'

workflow PREPARE_GENOME {

    take:
        fasta // channel (mandatory): [ val(meta2), path(index) ]

    main:
        ch_versions = Channel.empty()
        ch_index = Channel.empty()

        //TODO: Leverage more the igenomes machinery, avoid specifying parameters if not needed

        if (params.aligner == 'bwamem') { //if aligner is bwa-mem

            if (params.index_genome) {
                BWA_INDEX((fasta.map{ it -> [[id:it[0].baseName], it] }))
                ch_versions = ch_versions.mix(BWA_INDEX.out.versions.first())


                }

            ch_index = params.fasta ? params.bwa ? Channel.value(file(params.bwa)).map{ it -> [[id:it[0].baseName], it] } : BWA_INDEX.out.index : Channel.value(file(params.genomes[params.genome].bwa)).collect()

            } else { // if aligner is bwa-mem2

            if (params.index_genome) {
                BWAMEM2_INDEX((fasta.map{ it -> [[id:it[0].baseName], it] }))
                ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions.first())
                }
            ch_index  = params.fasta ? params.bwamem2  ? Channel.value(file(params.bwamem2))  : BWAMEM2_INDEX.out.index : Channel.value(file(params.genomes[params.genome].bwamem2)).collect()

            }

        if (params.index_genome) {
            SAMTOOLS_FAIDX(fasta.map{ it -> [[id:it[0].baseName], it] })
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())
        }
        ch_fasta_fai = params.fasta ? params.fasta_fai  ? Channel.value(file(params.fasta_fai))  : SAMTOOLS_FAIDX.out.fai : Channel.value(file(params.genomes[params.genome].fasta_fai)).collect()

    emit:
        index = ch_index.collect()
        fasta_fai = ch_fasta_fai.collect()
        versions = ch_versions
}
