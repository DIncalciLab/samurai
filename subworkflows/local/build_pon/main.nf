
include { HMMCOPY_READCOUNTER as HMMCOPY_READCOUNTER_PON   } from "../../../modules/nf-core//hmmcopy/readcounter/main"
include { ICHORCNA_GENERATE_PON                            } from "../../../modules/local/ichorcna/create_pon/main"
include { SAMBAMBA_FILTER                                  } from "../../../modules/local//sambamba/filterfragment/main"
include { WISECONDORX_CONVERT as NORMAL_CONVERT            } from '../../../modules/local/wisecondorx/convert/main'  
include { WISECONDORX_NEWREF                               } from '../../../modules/local/wisecondorx/newref/main'  

workflow BUILD_PON {
    take:
        normal_dir
        
    main:
        ch_versions = Channel.empty()

        bam_files = Channel
                        .fromFilePairs("${normal_dir}/*.{bam,bam.bai}",
                                checkIfExists: true)
                        .map {
                            meta, file ->
                            fmeta = [:]
                            fmeta.id = meta
                            tuple(fmeta, file[0], file[1])
                        }
        if (params.filter_bam_pon) {
                SAMBAMBA_FILTER(bam_files)
                bam_for_pon = SAMBAMBA_FILTER.out.filtered_bam
            } else {
                // Remove the BAM index for compatibility with the ReadCounter workflow
                bam_for_pon = bam_files
            }

        if (params.ichorcna) {
            
            HMMCOPY_READCOUNTER_PON(bam_for_pon)

            wigfiles = HMMCOPY_READCOUNTER_PON.out.wig.map {
                it ->
                it[1]
            }
            gc_wig              = Channel.value(params.gc_wig)
            map_wig             = Channel.value(params.map_wig)
            centromere          = Channel.value(params.centromere)
            reptime_file        = Channel.value(params.reptime_file)
            
            ICHORCNA_GENERATE_PON(wigfiles.collect(),
                                gc_wig, map_wig, centromere, reptime_file)

            normal_panel = ICHORCNA_GENERATE_PON.out.pon_file     
            ch_versions = ch_versions.mix(ICHORCNA_GENERATE_PON.out.versions)
                           
        }

        if (params.wisecondorx) {

                NORMAL_CONVERT(bam_files)
                normal_npz = NORMAL_CONVERT.out.npz.collect{ sample, npz_file -> file(npz_file) }
                ch_versions = ch_versions.mix(NORMAL_CONVERT.out.versions)
                WISECONDORX_NEWREF(normal_npz)
                normal_panel = WISECONDORX_NEWREF.out.npz_reference    
                ch_versions = ch_versions.mix(WISECONDORX_NEWREF.out.versions)
        }


    emit:
        normal_panel = normal_panel
        versions = ch_versions
}