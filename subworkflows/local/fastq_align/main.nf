//
// Alignment with BWA-MEM
//

// Include nf-core modules
include { BWA_MEM                     } from '../../../modules/nf-core/bwa/mem/main' 
include { BWAMEM2_MEM                 } from '../../../modules/nf-core/bwamem2/mem/main'




// Workflow

workflow FASTQ_ALIGN {

    take:
        ch_reads        // channel (mandatory): [ val(meta), [ path(reads) ] ]
        ch_index        // channel (mandatory): [ val(meta2), path(index) ]
        sort            // boolean (mandatory): true or false
        

    main: 
        ch_versions = Channel.empty()
        ch_bam = Channel.empty()

        if (params.aligner == 'bwa-mem') { //if aligner is bwa-mem
          
            BWA_MEM(ch_reads, ch_index, sort) 
            ch_bam = ch_bam.mix(BWA_MEM.out.bam)
            ch_versions = ch_versions.mix( BWA_MEM.out.versions.first() )

        } else { // if aligner is bwa-mem2

            BWAMEM2_MEM(ch_reads, ch_index, sort) 
            ch_bam = ch_bam.mix( BWAMEM2_MEM.out.bam )
            ch_versions = ch_versions.mix( BWAMEM2_MEM.out.versions.first() )
        }
               
    emit:
        bam = ch_bam                                // channel: [ val(meta), path(bam) ]   
        
        versions = ch_versions                          // channel: [ path(versions.yml) ]

}
        