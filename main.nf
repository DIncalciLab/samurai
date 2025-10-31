#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dincalcilab/samurai
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/dincalcilab/samurai
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMURAI                 } from './workflows/samurai'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_samurai_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_samurai_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_samurai_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// They need to be there before imports, cf. #19
params.fasta = getGenomeAttribute('fasta')
params.fai   = getGenomeAttribute('fasta_fai')
params.dict  = getGenomeAttribute('dict')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


//
// WORKFLOW: Run main dincalcilab/samurai analysis pipeline
//
workflow DINCALCILAB_SAMURAI {
    take:
    samplesheet
    aligner
    analysis_type
    genome
    fasta
    fai
    dict
    genome_index
    caller
    binsize
    pon_path
    options

    main:

    SAMURAI(
        samplesheet,
        aligner,
        analysis_type,
        genome,
        fasta,
        fai,
        dict,
        genome_index,
        caller,
        binsize,
        pon_path,
        options,
    )

    emit:
    multiqc_report = SAMURAI.out.multiqc_report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {


    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
    )

    def options = [
        build_pon: params.build_pon,
        index_genome: params.index_genome,
        run_fastp: params.run_fastp,
        compute_signatures: params.analysis_type != "align_only" ? params.compute_signatures : false,
        run_gistic: params.analysis_type != "align_only" ? params.run_gistic : false,
    ]

    genome = channel.value(params.genome)
    caller = params.analysis_type == "align_only" ? channel.value("none") : channel.value(params.caller)
    pon_path = params.pon_path && params.build_pon ? channel.value(params.pon_path) : []
    analysis_type = channel.value(params.analysis_type)
    binsize = channel.value(params.binsize)

    if (params.aligner == "bwamem") {
        real_aligner = "bwa"
    }
    else if (params.aligner == "bwamem2") {
        real_aligner = 'bwamem2'
    }
    else {
        real_aligner = ''
    }

    if (!params.aligner_index && !params.igenomes_ignore) {
        ch_index = [
            ["id": "aligner"],
            file(getGenomeAttribute(real_aligner)),
        ]
    }
    else {
        ch_index = [[], []]
    }

    // pass build_index as param here
    // pass size_selection as param here

    DINCALCILAB_SAMURAI(
        PIPELINE_INITIALISATION.out.samplesheet,
        real_aligner,
        analysis_type,
        genome,
        PIPELINE_INITIALISATION.out.fasta,
        PIPELINE_INITIALISATION.out.fai,
        PIPELINE_INITIALISATION.out.dict,
        ch_index,
        caller,
        binsize,
        pon_path,
        options,
    )

    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        DINCALCILAB_SAMURAI.out.multiqc_report,
    )
}
