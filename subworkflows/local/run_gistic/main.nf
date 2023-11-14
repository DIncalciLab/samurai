include { RUN_GISTIC2                                  } from "../../../modules/local/gistic2/main"

workflow RUN_GISTIC {
    take:
        caller
        gistic_file

    main:
        ch_versions = Channel.empty()
        
        switch(caller) {
            case "ichorcna":
                RUN_GISTIC2(gistic_file)                
                break
            case "wisecondorx":
                seg_file = file("${params.outdir}/wisecondorx/all_segments_wisecondorx_gistic.seg")
                RUN_GISTIC2(gistic_file)                
                break
            case "ascat_sc":
                seg_file = file("${params.outdir}/ascat_sc/all_segments_ascat_sc_gistic.seg")
                RUN_GISTIC2(gistic_file)                
                break    
            default:
                error "Unsupported GISTIC analysis for caller ${caller}"

        ch_versions = ch_versions.mix(RUN_GISTIC2.out.versions)
        }

    emit:
        versions = ch_versions
}

//TO DO: ADD POST PROCESSING OF GISTIC OUTPUT IN THE SUBWORKFLOW