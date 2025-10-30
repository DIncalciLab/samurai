include { HMMCOPY_READCOUNTER as HMMCOPY_READCOUNTER_PON } from "../../../modules/nf-core//hmmcopy/readcounter/main"
include { ICHORCNA_GENERATE_PON                          } from "../../../modules/local/ichorcna/create_pon/main"
include { SAMBAMBA_FILTER                                } from "../../../modules/local//sambamba/filterfragment/main"
include { WISECONDORX_CONVERT as NORMAL_CONVERT          } from '../../../modules/nf-core/wisecondorx/convert/main'
include { WISECONDORX_NEWREF                             } from '../../../modules/nf-core/wisecondorx/newref/main'

workflow BUILD_PON {
    take:
    normal_dir
    caller
    fasta
    fai

    main:

    ch_versions = channel.empty()
    ch_bam_files = channel.fromFilePairs(
            "${normal_dir}/*.bam{,.bai}",
            checkIfExists: true
        ).ifEmpty { error("No BAM or BAI files found at ${normal_dir}") }.map { meta, file ->
            def fmeta = [:]
            fmeta.id = meta
            tuple(fmeta, file[0], file[1])
        }

    if (caller == "ichorcna") {
        // FIXME: We shouldn't depend on parameters here
        if (params.filter_bam_pon) {
            SAMBAMBA_FILTER(ch_bam_files)
            ch_bam_for_pon = SAMBAMBA_FILTER.out.filtered_bam
        }
        else {
            // Remove the BAM index for compatibility with the ReadCounter workflow
            ch_bam_for_pon = ch_bam_files
        }

        HMMCOPY_READCOUNTER_PON(ch_bam_for_pon)
        wigfiles = HMMCOPY_READCOUNTER_PON.out.wig
            .map { _meta, wigfile ->
                wigfile
            }
            .collect()

        // FIXME: ship these with the pipeline
        // Use files and not values to avoid hangs (#20)
        gc_wig = file(params.ichorcna_gc_wig)
        map_wig = file(params.ichorcna_map_wig)
        reptime_file = file(params.ichorcna_reptime_wig)
        centromere = file(params.ichorcna_centromere_file)

        ICHORCNA_GENERATE_PON(
            wigfiles,
            gc_wig,
            map_wig,
            centromere,
            reptime_file,
        )

        normal_panel = ICHORCNA_GENERATE_PON.out.pon_file
        ch_versions = ch_versions.mix(ICHORCNA_GENERATE_PON.out.versions)
    }
    else if (caller == "wisecondorx") {
        NORMAL_CONVERT(ch_bam_files, fasta, fai)
        ch_versions = ch_versions.mix(NORMAL_CONVERT.out.versions)
        WISECONDORX_NEWREF(NORMAL_CONVERT.out.npz.map {meta, npz ->
            def new_meta = meta + [id: "joined"]
            [new_meta, npz]
            }.groupTuple())
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
