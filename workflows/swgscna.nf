/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'
include { fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

// Validate input parameters
WorkflowSwgscna.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def checkPathParamList = [ params.input,
                           params.multiqc_config,
                           params.fasta,
                           params.fai,
                           params.dict,
                           params.qdnaseq_bin_data
                           ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { INPUT_CHECK                } from '../subworkflows/local/input_check'
include { PREPARE_GENOME             } from '../subworkflows/local/prepare_genome/main'
include { SOLID_BIOPSY               } from '../subworkflows/local/solid_biopsy/main'
include { SIZE_SELECTION             } from '../subworkflows/local/size_selection/main'
include { LIQUID_BIOPSY              } from '../subworkflows/local/liquid_biopsy/main'
//include { COMPUTE_SIGNATURES         } from '../subworkflows/local/compute_signatures/main'
include { RUN_GISTIC                 } from '../subworkflows/local/run_gistic/main'

include { FASTA_INDEX_DNA            } from '../subworkflows/nf-core/fasta_index_dna/main'
include { FASTQ_ALIGN_DNA            } from '../subworkflows/nf-core/fastq_align_dna/main'
include { BAM_MARKDUPLICATES_PICARD  } from '../subworkflows/nf-core/bam_markduplicates_picard/main'
include { BAM_QC_PICARD              } from '../subworkflows/nf-core/bam_qc_picard/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { CIN_SIGNATURE_QUANTIFICATION} from '../modules/local/cin_signature_quantification/main'

include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

//
// SUBWORKFLOWS nf-core
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SWGSCNA {

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    binsize = Channel.value(params.binsize)
    genome = Channel.value(params.genome)
    caller = Channel.value(params.caller)

    ch_fasta = Channel.value(
        [["id": "fasta"], file(params.fasta, checkIfExists: true)]
    )
    ch_fai = Channel.value(
        [["id": "fai"], file(params.fai, checkIfExists: true)]
    )

    ch_dict = Channel.value(
        [["id": "dict"], file(params.dict, checkIfExists: true)]
    )

    // Dummy channel for indexing
    ch_liftover = Channel.value(
        [["id": "liftover"], []]
    )
    //FIXME: Differentiate between BAM / FASTQ
    ch_input = Channel.fromSamplesheet("input").map {
        meta, fastq1, fastq2 ->
            if(fastq2) {
                [meta, [fastq1, fastq2]]
            } else {
                [meta, [fastq1]]
            }
    }

    if (params.aligner) {
        FASTQC ( ch_input )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())

        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

        // SUBWORKFLOW: FASTQ_ALIGN_DNA

        if (params.index_genome) {
            FASTA_INDEX_DNA(ch_fasta, ch_liftover, params.aligner)
            ch_index = FASTA_INDEX_DNA.out.index
            ch_versions = ch_versions.mix(FASTA_INDEX_DNA.out.versions.first())
        } else {
            if (!params.aligner_index && !params.igenomes_ignore) {
                ch_index = [["id": "aligner"], file(WorkflowMain.getGenomeAttribute(params, params.aligner))]
            } else {
                ch_index = [["id": "aligner"], file(params.aligner_index)]
            }
        }

        FASTQ_ALIGN_DNA (
            ch_input,
            ch_index,
            params.aligner,
            true // sort_bam
        )
        ch_versions = ch_versions.mix(FASTQ_ALIGN_DNA.out.versions.first())

        // Mark duplicates after alignment
        BAM_MARKDUPLICATES_PICARD (
            FASTQ_ALIGN_DNA.out.bam,
            ch_fasta,
            ch_fai
        )
        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions.first())

        ch_multiqc_files = ch_multiqc_files.mix(
            BAM_MARKDUPLICATES_PICARD.out.metrics.collect{
                meta, metrics -> metrics
            }
        )
        ch_multiqc_files = ch_multiqc_files.mix(
            BAM_MARKDUPLICATES_PICARD.out.flagstat.collect{
                meta, metrics -> metrics
            }
        )
        ch_multiqc_files = ch_multiqc_files.mix(
                BAM_MARKDUPLICATES_PICARD.out.idxstats.collect{
                    meta, metrics -> metrics
                }
        )

        ch_bam_bai = BAM_MARKDUPLICATES_PICARD.out.bam
            .join(BAM_MARKDUPLICATES_PICARD.out.bai, by: [0], remainder: true)
            .join(BAM_MARKDUPLICATES_PICARD.out.csi, by: [0], remainder: true)
            .map {
                meta, bam, bai, csi ->
                    if (bai) {
                        [ meta, bam, bai ]
                    } else {
                        [ meta, bam, csi ]
                }
        }

        // QC metrics about alignment (coverage, etc.)
        BAM_QC_PICARD(
            // sWGS has neither baits nor targets, 3rd and 4th args are thus empty
            ch_bam_bai.map {
                meta, bam, bai -> [meta, bam, bai, [], []]
            },
            ch_fasta,
            ch_fai,
            ch_dict
        )

        ch_versions = ch_versions.mix(
            BAM_QC_PICARD.out.versions.first()
        )
        ch_multiqc_files = ch_multiqc_files.mix(
                BAM_QC_PICARD.out.coverage_metrics.collect{
                meta, metrics -> metrics
            }
        )
        ch_multiqc_files = ch_multiqc_files.mix(
            BAM_QC_PICARD.out.multiple_metrics.collect{
                meta, metrics -> metrics
            }
        )

    } else {
        ch_bam_bai = INPUT_CHECK.out.reads //TODO: switch to nf-validation
    }


    // CN Calling

    switch(params.analysis_type) {
        case "solid_biopsy":
            SOLID_BIOPSY(ch_bam_bai, params.caller, binsize, Channel.value(params.genome))
            ch_versions = ch_versions.mix(SOLID_BIOPSY.out.versions.first())
            ch_multiqc_files = ch_multiqc_files.mix(SOLID_BIOPSY.out.summary.collect())
            break
        case "liquid_biopsy":
            if (params.size_selection && params.caller == "ichorcna") {
                SIZE_SELECTION(ch_bam_bai, file(params.fasta))
                ch_versions = ch_versions.mix(SIZE_SELECTION.out.versions.first())

                ch_multiqc_files = ch_multiqc_files.mix(SIZE_SELECTION.out.stats_pre.map{
                        meta, file -> return file })
                ch_multiqc_files = ch_multiqc_files.mix(SIZE_SELECTION.out.stats_post.map{
                        meta, file -> return file })

                ch_analysis = SIZE_SELECTION.out.bam.join(SIZE_SELECTION.out.bai, by: [0], remainder: true)
                            .map {
                                meta, bam, bai -> [ meta, bam, bai ]
                            }
                } else {
                    ch_analysis = ch_bam_bai
                }
            LIQUID_BIOPSY(ch_analysis, params.caller)
            ch_versions = ch_versions.mix(LIQUID_BIOPSY.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(LIQUID_BIOPSY.out.summary.collect())
            break
        default:
            error "Uknown / unsupported analysis ${analysis_type}"
    }

    // Compute CN Signatures if specified, default: false
    if (params.compute_signatures) {
        CIN_SIGNATURE_QUANTIFICATION(LIQUID_BIOPSY.out.signature_file)
        ch_versions = ch_versions.mix(CIN_SIGNATURE_QUANTIFICATION.out.versions)
    }

    // Run GISTIC if specified, default: false
    if (params.run_gistic) {
        RUN_GISTIC(params.caller, LIQUID_BIOPSY.out.gistic_file)
        ch_versions = ch_versions.mix(RUN_GISTIC.out.versions)
    }

    // Software versions
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowSwgscna.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowSwgscna.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
