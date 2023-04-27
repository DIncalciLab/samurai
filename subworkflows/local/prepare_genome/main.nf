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

        if (params.index_genome) {

            if (params.aligner == 'bwa-mem') { //if aligner is bwa-mem

                BWA_INDEX((fasta.map{ it -> [[id:it[0].baseName], it] }))
                ch_index = ch_index.mix(BWA_INDEX.out.index)
                ch_versions = ch_versions.mix(BWA_INDEX.out.versions.first())

            } else { // if aligner is bwa-mem2 

                BWAMEM2_INDEX((fasta.map{ it -> [[id:it[0].baseName], it] }))
                ch_index = ch_index.mix(BWAMEM2_INDEX.out.index)
                ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions.first())
            }
        } else {
            ch_index = ch_index.mix(ch_fasta)
        }

                
        SAMTOOLS_FAIDX(fasta.map{ it -> [[id:it[0].baseName], it] })
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())
        
    emit:
        index = ch_index.map{ meta, index -> [index] }.collect()     
        fasta_fai = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }
        //gzi = SAMTOOLS_FAIDX.out.gzi
        versions = ch_versions     
}