
include { HMMCOPY_READCOUNTER as HMMCOPY_READCOUNTER_PON   } from "../../../modules/nf-core//hmmcopy/readcounter/main"
include { ICHORCNA_GENERATE_PON                            } from "../../../modules/local/ichorcna/create_pon/main"
include { SAMBAMBA_FILTER                                  } from "../../../modules/local//sambamba/filterfragment/main"
include { WISECONDORX_CONVERT as NORMAL_CONVERT            } from '../../../modules/local/wisecondorx/convert/main'
include { WISECONDORX_NEWREF                               } from '../../../modules/local/wisecondorx/newref/main'

workflow BUILD_PON {
    take:
        normal_dir
        caller

    main:
        ch_versions = Channel.empty()
        ch_bam_files = Channel
                        .fromFilePairs("${normal_dir}/*.{bam,bam.bai}",
                                checkIfExists: true)
                        .map {
                            meta, file ->
                            fmeta = [:]
                            fmeta.id = meta
                            tuple(fmeta, file[0], file[1])
                        }

        switch(caller) {
            case "ichorcna":

                 // FIXME: We shouldn't depend on parameters here
                if (params.filter_bam_pon) {
                    SAMBAMBA_FILTER(ch_bam_files)
                    ch_bam_for_pon = SAMBAMBA_FILTER.out.filtered_bam
                } else {
                // Remove the BAM index for compatibility with the ReadCounter workflow
                ch_bam_for_pon = ch_bam_files
                }

                HMMCOPY_READCOUNTER_PON(ch_bam_for_pon)
                wigfiles = HMMCOPY_READCOUNTER_PON.out.wig.map {
                    it ->
                    it[1]
                }

                gc_wig              = Channel.value(params.ichorcna_gc_wig)
                map_wig             = Channel.value(params.ichorcna_map_wig)
                reptime_file        = Channel.value(params.ichorcna_reptime_wig)
                centromere          = Channel.value(params.ichorcna_centromere_file)

                ICHORCNA_GENERATE_PON(wigfiles.collect(),
                                    gc_wig, map_wig, centromere, reptime_file)

                normal_panel = ICHORCNA_GENERATE_PON.out.pon_file
                ch_versions = ch_versions.mix(ICHORCNA_GENERATE_PON.out.versions)
                break
            case "wisecondorx":
                NORMAL_CONVERT(ch_bam_for_pon)
                ch_normal_npz = NORMAL_CONVERT.out.npz.collect{ sample, npz_file -> file(npz_file) }
                ch_versions = ch_versions.mix(NORMAL_CONVERT.out.versions)
                WISECONDORX_NEWREF(ch_normal_npz)
                normal_panel = WISECONDORX_NEWREF.out.npz_reference
                ch_versions = ch_versions.mix(WISECONDORX_NEWREF.out.versions)
                break
            default:
                error "Unknown/unsupported caller ${caller}"
        }

    emit:
        normal_panel = normal_panel
        versions = ch_versions
}
