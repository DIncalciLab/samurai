include { HMMCOPY_READCOUNTER as HMMCOPY_READCOUNTER_PON } from "../../../modules/nf-core//hmmcopy/readcounter/main"
include { ICHORCNA_CREATEPON                             } from '../../../modules/nf-core/ichorcna/createpon/main'
include { SAMBAMBA_FILTER                                } from "../../../modules/local//sambamba/filterfragment/main"
include { WISECONDORX_CONVERT as NORMAL_CONVERT          } from '../../../modules/nf-core/wisecondorx/convert/main'
include { WISECONDORX_NEWREF                             } from '../../../modules/nf-core/wisecondorx/newref/main'

workflow BUILD_PON {
    take:
    normal_dir
    caller
    fasta
    fai
    gc_wig
    map_wig
    reptime
    centromere
    filter_bam_pon

    main:

    ch_versions = channel.empty()
    ch_bam_files = channel.fromFilePairs(
            "${normal_dir}/*.bam{,.bai}",
            checkIfExists: true
        )
        .ifEmpty { error("No BAM or BAI files found at ${normal_dir}") }
        .map { meta, file ->
            def fmeta = [:]
            fmeta.id = meta
            tuple(fmeta, file[0], file[1])
        }

    if (caller == "ichorcna") {
        // FIXME: We shouldn't depend on parameters here
        if (filter_bam_pon) {
            SAMBAMBA_FILTER(ch_bam_files)
            ch_bam_for_pon = SAMBAMBA_FILTER.out.filtered_bam
            ch_versions = ch_versions.mix(SAMBAMBA_FILTER.out.versions)
        }
        else {
            // Remove the BAM index for compatibility with the ReadCounter workflow
            ch_bam_for_pon = ch_bam_files
        }

        HMMCOPY_READCOUNTER_PON(ch_bam_for_pon, fasta)
        wigfiles = HMMCOPY_READCOUNTER_PON.out.wig
            .map { _meta, wigfile ->
                wigfile
            }
            .collect()

        ch_versions = ch_versions.mix(HMMCOPY_READCOUNTER_PON.out.versions)

        ICHORCNA_CREATEPON(
            wigfiles,
            gc_wig,
            map_wig,
            centromere,
            reptime,
            [], /* exons */
        )

        normal_panel = ICHORCNA_CREATEPON.out.rds
        ch_versions = ch_versions.mix(ICHORCNA_CREATEPON.out.versions)
    }
    else if (caller == "wisecondorx") {
        NORMAL_CONVERT(ch_bam_files, fasta, fai)
        ch_versions = ch_versions.mix(NORMAL_CONVERT.out.versions)
        WISECONDORX_NEWREF(
            NORMAL_CONVERT.out.npz.map { meta, npz ->
                def new_meta = meta + [id: "joined"]
                [new_meta, npz]
            }.groupTuple()
        )
        normal_panel = WISECONDORX_NEWREF.out.npz
        ch_versions = ch_versions.mix(WISECONDORX_NEWREF.out.versions)
    }
    else {
        error("Unknown/unsupported caller ${caller}")
    }

    emit:
    normal_panel = normal_panel
    versions     = ch_versions
}
