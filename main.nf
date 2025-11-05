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
    run_fastp
    build_pon
    normal_panel
    index_genome
    run_gistic
    size_selection
    wisecondor_blacklist
    ichorcna_gc_wig
    ichorcna_map_wig
    ichorcna_centromere
    ichorcna_reptime_file
    ascat_sc_predict_refit

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
        run_fastp,
        build_pon,
        normal_panel,
        index_genome,
        run_gistic,
        size_selection,
        wisecondor_blacklist,
        ichorcna_gc_wig,
        ichorcna_map_wig,
        ichorcna_centromere,
        ichorcna_reptime_file,
        ascat_sc_predict_refit,
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

    if (!params.fasta) {
        error("Error: a reference FASTA file was not provided.")
    }

    if (!params.fai) {
        error("Error: a FASTA index was not provided.")
    }

    if (!params.dict) {
        error("Error: a FASTA sequence dictionary was not provided.")
    }

    fasta = params.fasta ? channel.fromPath(params.fasta).map { it -> [[id: it.baseName], it] }.collect() : channel.empty()
    fai = params.fai ? channel.fromPath(params.fai).map { it -> [[id: it.baseName], it] }.collect() : channel.empty()
    dict = params.dict ? channel.fromPath(params.fai).map { it -> [[id: it.baseName], it] }.collect() : channel.empty()

    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
    )

    ch_samplesheet = PIPELINE_INITIALISATION.out.samplesheet.map { it ->
        def meta = it[0]
        def read_group = "\"@RG\\tID:${meta.sample}\\tPU:1\\tSM:${meta.sample}\\tLB:${meta.sample}\\tDS:${params.fasta}\\tPL:Illumina\""
        meta += [read_group: read_group]
        it[0] = meta
        it
    }

    run_gistic = params.analysis_type != "align_only" ? params.run_gistic : false

    genome = channel.value(params.genome)
    caller = params.analysis_type == "align_only" ? "none" : params.caller
    pon_path = params.pon_path && params.build_pon ? channel.value(params.pon_path) : channel.empty()
    analysis_type = params.analysis_type
    binsize = params.binsize
    normal_panel = params.normal_panel ? channel.fromPath(params.normal_panel, checkIfExists: true) : channel.empty()

    if (params.aligner == "bwamem") {
        real_aligner = "bwa"
    }
    else if (params.aligner == "bwamem2") {
        real_aligner = 'bwamem2'
    }
    else {
        real_aligner = ''
    }

    if (!params.aligner_index && !params.igenomes_ignore && real_aligner) {
        ch_index = [
            ["id": "aligner"],
            file(getGenomeAttribute(real_aligner)),
        ]
    }
    else {
        ch_index = [[], []]
    }

    ichorcna_gc_wig = params.ichorcna_gc_wig ? file(params.ichorcna_gc_wig, checkIfExists: true) : []
    ichorcna_map_wig = params.ichorcna_map_wig ? file(params.ichorcna_map_wig, checkIfExists: true) : []
    ichorcna_centromere = params.ichorcna_centromere_file ? file(params.ichorcna_centromere_file, checkIfExists: true) : []
    ichorcna_reptime_file = params.ichorcna_reptime_wig ? file(params.ichorcna_reptime_wig, checkIfExists: true) : []

    if (caller == "ichorcna" && (!ichorcna_gc_wig || !ichorcna_map_wig)) {
        error("ichorCNA calling requires a GC WIG and a mappability WIG")
    }

    wisecondor_blacklist = params.wisecondorx_blacklist ? channel.fromPath(params.wisecondorx_blacklist, checkIfExists: true).map { blacklist -> [[id: "blacklist"], blacklist] } : [[], []]

    // pass build_index as param here
    // pass size_selection as param here

    DINCALCILAB_SAMURAI(
        ch_samplesheet,
        real_aligner,
        analysis_type,
        genome,
        fasta,
        fai,
        dict,
        ch_index,
        caller,
        binsize,
        pon_path,
        params.run_fastp,
        params.build_pon,
        normal_panel,
        params.index_genome,
        run_gistic,
        params.size_selection,
        wisecondor_blacklist,
        ichorcna_gc_wig,
        ichorcna_map_wig,
        ichorcna_centromere,
        ichorcna_reptime_file,
        params.ascat_sc_predict_refit,
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
