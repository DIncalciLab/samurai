/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap             } from 'plugin/nf-schema'
include { paramsSummaryMultiqc         } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML       } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText       } from '../subworkflows/local/utils_samurai_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { CIN_SIGNATURE_QUANTIFICATION } from '../modules/local/cin_signature_quantification/main'
include { SOLID_BIOPSY                 } from '../subworkflows/local/solid_biopsy/main'
include { SIZE_SELECTION               } from '../subworkflows/local/size_selection/main'
include { LIQUID_BIOPSY                } from '../subworkflows/local/liquid_biopsy/main'
include { RUN_GISTIC                   } from '../subworkflows/local/run_gistic/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { FASTQC                       } from '../modules/nf-core/fastqc/main'
include { MULTIQC                      } from '../modules/nf-core/multiqc/main'
include { SAMTOOLS_INDEX               } from '../modules/nf-core/samtools/index/main'

//
// SUBWORKFLOWS nf-core
//
include { FASTQ_TRIM_FASTP_FASTQC      } from '../subworkflows/nf-core/fastq_trim_fastp_fastqc/main'
include { FASTA_INDEX_DNA              } from '../subworkflows/nf-core/fasta_index_dna/main'
include { FASTQ_ALIGN_DNA              } from '../subworkflows/nf-core/fastq_align_dna/main'
include { BAM_MARKDUPLICATES_PICARD    } from '../subworkflows/nf-core/bam_markduplicates_picard/main'
include { BAM_QC_PICARD                } from '../subworkflows/nf-core/bam_qc_picard/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SAMURAI {
    take:
    ch_input // channel: samplesheet read in from --input
    aligner // value: supplied aligner
    analysis_type
    genome
    ch_fasta
    ch_fai
    ch_dict
    ch_index
    caller
    binsize
    ch_pon_path
    run_fastp
    build_pon
    ch_normal_panel
    index_genome
    run_gistic
    size_selection
    ch_blacklist
    ch_gc_wig
    ch_map_wig
    ch_centromere
    ch_reptime
    ascat_predict_refit

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    // Dummy channel for indexing
    ch_liftover = channel.value(
        [["id": "liftover"], []]
    )

    if (aligner) {
        skip_fastp = run_fastp ? false : true

        ch_fastq = ch_input.map { meta, fastq -> [meta, fastq, []] }

        FASTQ_TRIM_FASTP_FASTQC(
            ch_fastq,
            false,
            false,
            false,
            skip_fastp,
            false,
        )
        ch_versions = ch_versions.mix(FASTQ_TRIM_FASTP_FASTQC.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(
            FASTQ_TRIM_FASTP_FASTQC.out.fastqc_raw_zip.collect { it -> it[1] }.ifEmpty([])
        )
        ch_multiqc_files = ch_multiqc_files.mix(
            FASTQ_TRIM_FASTP_FASTQC.out.trim_json.collect { it -> it[1] }.ifEmpty([])
        )
        ch_multiqc_files = ch_multiqc_files.mix(
            FASTQ_TRIM_FASTP_FASTQC.out.fastqc_trim_zip.collect { it -> it[1] }.ifEmpty([])
        )

        ch_input = FASTQ_TRIM_FASTP_FASTQC.out.reads

        // SUBWORKFLOW: FASTQ_ALIGN_DNA

        if (index_genome) {
            FASTA_INDEX_DNA(ch_fasta, ch_liftover, aligner)
            ch_index = FASTA_INDEX_DNA.out.index
            ch_versions = ch_versions.mix(FASTA_INDEX_DNA.out.versions.first())
        }

        FASTQ_ALIGN_DNA(
            ch_input,
            ch_index,
            ch_fasta,
            aligner,
            true,
        )
        ch_versions = ch_versions.mix(FASTQ_ALIGN_DNA.out.versions.first())

        // Mark duplicates after alignment
        BAM_MARKDUPLICATES_PICARD(
            FASTQ_ALIGN_DNA.out.bam,
            ch_fasta,
            ch_fai,
        )
        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions.first())

        ch_multiqc_files = ch_multiqc_files.mix(
            BAM_MARKDUPLICATES_PICARD.out.metrics.collect { _meta, metrics ->
                metrics
            }
        )
        ch_multiqc_files = ch_multiqc_files.mix(
            BAM_MARKDUPLICATES_PICARD.out.flagstat.collect { _meta, metrics ->
                metrics
            }
        )
        ch_multiqc_files = ch_multiqc_files.mix(
            BAM_MARKDUPLICATES_PICARD.out.idxstats.collect { _meta, metrics ->
                metrics
            }
        )

        ch_bam_bai = BAM_MARKDUPLICATES_PICARD.out.bam
            .join(BAM_MARKDUPLICATES_PICARD.out.bai, by: [0], remainder: true)
            .join(BAM_MARKDUPLICATES_PICARD.out.csi, by: [0], remainder: true)
            .map { meta, bam, bai, csi ->
                if (bai) {
                    [meta, bam, bai]
                }
                else {
                    [meta, bam, csi]
                }
            }
    }
    else {
        // Just index the files and we're good to go
        SAMTOOLS_INDEX(ch_input)
        ch_bam_bai = ch_input
            .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
            .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
            .map { meta, bam, bai, csi ->
                if (bai) {
                    [meta, bam, bai]
                }
                else {
                    [meta, bam, csi]
                }
            }
    }

    // QC metrics about alignment (coverage, etc.)
    BAM_QC_PICARD(
        ch_bam_bai.map { meta, bam, bai ->
            [meta, bam, bai, [], []]
        },
        ch_fasta,
        ch_fai,
        ch_dict,
        [[], []] /* fasta_gzi */
    )

    ch_versions = ch_versions.mix(
        BAM_QC_PICARD.out.versions.first()
    )
    ch_multiqc_files = ch_multiqc_files.mix(
        BAM_QC_PICARD.out.coverage_metrics.collect { _meta, metrics ->
            metrics
        }
    )
    ch_multiqc_files = ch_multiqc_files.mix(
        BAM_QC_PICARD.out.multiple_metrics.collect { _meta, metrics ->
            metrics
        }
    )

    // CN Calling

    ch_bam_type = ch_bam_bai.branch{ meta, _bam, _bai ->
        normal: meta.status == "normal"
        tumor: true // fallback if not specified
    }

    ch_ponfiles = build_pon ? ch_bam_type.normal.ifEmpty([[], []]) : [[], []]

    if (analysis_type == "solid_biopsy") {
        SOLID_BIOPSY(
            ch_bam_type.tumor,
            caller,
            binsize,
            genome,
            ascat_predict_refit,
            build_pon,
            ch_ponfiles,
            ch_normal_panel,
            ch_gc_wig,
            ch_map_wig,
            ch_centromere,
            ch_reptime,
            ch_fasta
        )
        gistic_file = SOLID_BIOPSY.out.gistic_file
        ch_versions = ch_versions.mix(SOLID_BIOPSY.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(SOLID_BIOPSY.out.summary.collect())
    }
    else if (analysis_type == "liquid_biopsy") {
        if (size_selection) {
            SIZE_SELECTION(ch_bam_type.tumor, ch_fasta)
            ch_versions = ch_versions.mix(SIZE_SELECTION.out.versions.first())

            ch_multiqc_files = ch_multiqc_files.mix(
                SIZE_SELECTION.out.stats_pre.map { _meta, file ->
                    return file
                }
            )
            ch_multiqc_files = ch_multiqc_files.mix(
                SIZE_SELECTION.out.stats_post.map { _meta, file ->
                    return file
                }
            )

            ch_analysis = SIZE_SELECTION.out.bam
                .join(SIZE_SELECTION.out.bai, by: [0], remainder: true)
                .map { meta, bam, bai ->
                    [meta, bam, bai]
                }
        }
        else {
            ch_analysis = ch_bam_type.tumor
        }

        LIQUID_BIOPSY(
            ch_analysis,
            caller,
            ch_fasta,
            ch_fai,
            ch_normal_panel,
            ch_gc_wig,
            ch_map_wig,
            ch_centromere,
            ch_reptime,
            build_pon,
            ch_ponfiles,
            ch_blacklist)
        gistic_file = LIQUID_BIOPSY.out.corrected_gistic_file
        ch_versions = ch_versions.mix(LIQUID_BIOPSY.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(LIQUID_BIOPSY.out.summary.collect())
    }
    else if (analysis_type != "align_only") {
        error("Unknown / unsupported analysis ${analysis_type}")
    }

    // Run GISTIC if specified
    if (run_gistic) {
        RUN_GISTIC(gistic_file, genome)
        ch_versions = ch_versions.mix(RUN_GISTIC.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(RUN_GISTIC.out.gistic_lesions)
        ch_multiqc_files = ch_multiqc_files.mix(RUN_GISTIC.out.gistic_broad_lesions)
        ch_multiqc_files = ch_multiqc_files.mix(RUN_GISTIC.out.chrom_plot)
    }

    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf-hrsc_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? channel.fromPath(params.multiqc_config, checkIfExists: true)
        : channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : channel.empty()

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
