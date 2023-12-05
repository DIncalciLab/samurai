
include { HMMCOPY_READCOUNTER as HMMCOPY_READCOUNTER_PON   } from "../../../modules/nf-core//hmmcopy/readcounter/main"
include { ICHORCNA_GENERATE_PON                            } from "../../../modules/local/ichorcna/create_pon/main"
include { SAMBAMBA_FILTER                                  } from "../../../modules/local//sambamba/filterfragment/main"

workflow ICHORCNA_PON {
    take:
        normal_dir
        gc_wig
        map_wig
        centromere
        reptime_file

    main:

        bam_files = Channel
                    .fromFilePairs("${normal_dir}/*.{bam,bam.bai}",
                            checkIfExists: true)
                    .map {
                        meta, file ->
                        fmeta = [:]
                        fmeta.id = meta
                        tuple(fmeta, file[0], file[1])
                    }

        // FIXME: We shouldn't depend on parameters here
        if (params.filter_bam_pon) {
            SAMBAMBA_FILTER(bam_files)
            bam_for_pon = SAMBAMBA_FILTER.out.filtered_bam
        } else {
            // Remove the BAM index for compatibility with the ReadCounter workflow
            bam_for_pon = bam_files
        }

        HMMCOPY_READCOUNTER_PON(bam_for_pon)

        wigfiles = HMMCOPY_READCOUNTER_PON.out.wig.map {
            it ->
            it[1]
        }

        ICHORCNA_GENERATE_PON(wigfiles.collect(),
                            gc_wig, map_wig, centromere, reptime_file)


    emit:
        normal_panel = ICHORCNA_GENERATE_PON.out.pon_file
        versions = ICHORCNA_GENERATE_PON.out.versions
}
