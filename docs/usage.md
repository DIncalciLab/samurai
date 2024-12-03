# Introduction

`dincalcilab/samurai` is a bioinformatics best-practice analysis pipeline for the analysis of shallow whole genome sequencing (sWGS) data for the identification of copy number alterations (CNAs). It supports a number of workflows depending on the nature of the samples (coming from tissues or other biological fluids like plasma). While SAMURAI was developed with cancer studies in mind, it is applicable to any field where DNA alterations need to be studied.

# Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with different columns depending on the sample you would like to analyze.

If you want to start your analysis from FASTQ files of the sample, you need to set up your `samplesheet.csv` with two **required** columns at least:

```console
sample,fastq_1
SAMPLE1,path/to/SAMPLE1.fastq.gz
SAMPLE2,path/to/SAMPLE2.fastq.gz
```

If you have paired-end sequencing samples, your `samplesheet.csv` will look like:

```console
sample,fastq_1,fastq_2
SAMPLE1,path/to/SAMPLE1_R1.fastq.gz,path/to/SAMPLE1_R2.fastq.gz
SAMPLE2,path/to/SAMPLE2_R1.fastq.gz,path/to/SAMPLE2_R2.fastq.gz
```

You may want to skip the alignment step and start from pre-aligned BAM files. In this case your `samplesheet.csv` will look like:

```console
sample,bam,gender
SAMPLE1,path/to/SAMPLE1.bam,female
SAMPLE2,path/to/SAMPLE2.bam,male
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1.                                                                                                                                    |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2.                                                                                                                                    |
| `bam`     | Full path to BAM file. Note: `bam` is _mutually exclusive_ with `fastq_1` or `fastq_2`.                                                                                                |

| Optional | Description                                                                                                  |
| -------- | ------------------------------------------------------------------------------------------------------------ |
| `gender` | Gender of the patient that may be used in the copy number analysis; this could be either `female` or `male`. |

> **NB:** You can either use FASTQ files, or BAM files in a samplesheet. Mixing FASTQ fles and BAM files in the same samplesheet **is not supported**.

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

# Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run dincalcilab/samurai --input samplesheet.csv --outdir <OUTDIR> --genome hg38 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

#### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull dincalcilab/samurai
```

#### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [dincalcilab/samurai releases page](https://github.com/dincalcilab/samurai/releases) and find the latest pipeline version - numeric only (eg. `1.0.4`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.0.4`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

# Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

#### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test_ichorcna,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is **_not_** recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test_ascat_sc`
  - A profile with a complete configuration for automated testing of `solid_biopsy` workflow with `ASCAT.sc`
  - Includes links to test data so needs no other parameters
- `test_ichorcna`
  - A profile with a complete configuration for automated
    testing of `liquid_biopsy` workflow with `ichorCNA`
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda **as a last resort** i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.

#### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

#### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

# Custom configuration

> **NB:** These options are specific of SAMURAI and use a _double_ hyphen (ex `--caller`).

## General SAMURAI options

| Parameter | Description                                                                           |
| --------- | ------------------------------------------------------------------------------------- |
| `input`   | Path to `samplesheet.csv` containing information about the samples in the experiment. |
| `outdir`  | The output directory where the results will be saved.                                 |

## Reference genome options

| Parameter | Description                                                                                                                                                                                                                                       |
| --------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `genome`  | If using a reference genome configured in the pipeline using iGenomes, use this parameter to give the ID for the reference. See the [nf-core website docs](https://nf-co.re/usage/reference_genomes) for more details. Default is `--genome hg38` |
| `fasta`   | Path to FASTA genome file. If you are using iGenomes resources, it is automatically retrieved.                                                                                                                                                    |
| `fai`     | Path to FASTA FAI index. If you are using iGenomes resources, it is automatically retrieved.                                                                                                                                                      |
| `dict`    | Path to the FASTA sequence dictionary. If you are using iGenomes resources, it is automatically retrieved.                                                                                                                                        |

## Alignment options

| Parameter       | Description                                                                                                                                                             | Possible Values              |
| --------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------- |
| `aligner`       | This parameter specifies which aligner to use for the alignment step.If you start from BAM files and you want to skip the alingment step, you can set `--aligner false` | `bwamem`, `bwamem2`, `false` |
| `index_genome`  | Wether to index the reference genome or not.                                                                                                                            | `true`, `false`              |
| `aligner_index` | Path to the aligner index basename (must be compatible).                                                                                                                |

## FASTQ Trimming specific options

| Parameter Name          | Description                                                    | Possible Values              |
| ----------------------- | -------------------------------------------------------------- | ---------------------------- |
| `run_fastp`             | Run fastp for quality and UMI trimming                         | `true`, `false`              |
| `fastp_min_read_length` | Minimum length for a read to be kept after trimming            |                              |
| `fastp_cut_window_size` | Size of the sliding window used for quality trimming           |                              |
| `fastp_min_quality`     | Minimum Phred score to keep a base                             |                              |
| `fastp_umi_loc`         | Location of UMIs in the read                                   | `read1`, `read2`, `per_read` |
| `fastp_umi_skip`        | How many bases to skip after a UMI                             |                              |
| `fastp_umi_length`      | Length of the UMI to trim                                      |                              |
| `fastp_max_trimmed_pct` | Maximum percentage of trimmed bases before discarding a read   |                              |
| `fastp_trim_poly_x`     | Trim poly-X stretches in a read                                |                              |
| `fastp_trim_poly_g`     | Trim poly-G stretches in a read (NextSeq sequencing artifacts) |                              |

## Copy Number calling options

### Common Options

| Parameter Name  | Description                     | Default Value  | Possible Values                                  |
| --------------- | ------------------------------- | -------------- | ------------------------------------------------ |
| `binsize`       | Size of the genomic bins in kbp | 30             | Any integer (bin size in kilobases)              |
| `caller`        | Caller for sWGS                 | `qdnaseq`      | `ichorcna`, `wisecondorx`, `ascat_sc`, `qdnaseq` |
| `analysis_type` | Subworkflow to use              | `solid_biopsy` | `liquid_biopsy`, `solid_biopsy`, `"align_only"`  |

### Solid Biopsy Options

#### _QDNAseq specific options_

| Parameter Name        | Description                                                  | Default Value | Possible Values     |
| --------------------- | ------------------------------------------------------------ | ------------- | ------------------- |
| `qdnaseq_bin_data`    | Optional parameter to an RDS file containing bin annotations |               | Any valid file path |
| `qdnaseq_paired_ends` | Whether reads are paired or not                              | `true`        | `true`, `false`     |

#### _ASCAT.sc specific options_

| Parameter Name                | Description                                               | Default Value | Possible Values            |
| ----------------------------- | --------------------------------------------------------- | ------------- | -------------------------- |
| `ascat_sc_predict_refit`      | Perform refitting to select the optimal solution          | `"TRUE"`      | `"TRUE"`, `"FALSE"`        |
| `ascat_sc_segmentation_alpha` | Alpha value to use for segmentation                       | 0.01          | Any positive number        |
| `ascat_sc_min_purity`         | Minimum purity to use for prediction, in fractional units | 0.01          | Any number between 0 and 1 |
| `ascat_sc_max_purity`         | Maximum purity to use for prediction, in fractional units | 1             | Any number between 0 and 1 |
| `ascat_sc_min_ploidy`         | Minimum ploidy to consider                                | 1.7           | Any number greater than 1  |
| `ascat_sc_max_ploidy`         | Maximum ploidy to consider                                | 5             | Any number greater than 1  |
| `ascat_sc_max_tumor_ploidy`   | Maximum tumor ploidy to consider                          | 5             | Any integer greater than 1 |

### Liquid Biopsy Options

#### _Common Liquid biopsy options_

| Parameter Name               | Description                                                              | Default Value | Possible Values      |
| ---------------------------- | ------------------------------------------------------------------------ | ------------- | -------------------- |
| `normal_panel`               | Path to the panel of normals to be used                                  |               |                      |
| `pon_path`                   | Path to BAM files to be used to build the panel of normals               |               |                      |
| `pon_name`                   | Name of the panel of normals to build                                    | `"PoN"`       |                      |
| `build_pon`                  | Whether to build a panel of normals or not                               | `false`       | `true`, `false`      |
| `selection_maxsize`          | Maximum insert size in bp to be included for size selection              | `150`         | Any positive integer |
| `plot_fragment_distribution` | Whether to plot fragment size distributions during size selection (slow) | `false`       | `true`, `false`      |

#### _ichorCNA specific options_

| Parameter Name                      | Description                                                                 | Default Value                         | Possible Values                                  |
| ----------------------------------- | --------------------------------------------------------------------------- | ------------------------------------- | ------------------------------------------------ |
| `ichorcna_genome_style`             | Genome style to be used by ichorCNA (NCBI or UCSC)                          | `UCSC`                                | `NCBI`, `UCSC`                                   |
| `ichorcna_readcounter_chrs`         | Chromosomes to be used by ichorCNA                                          | `chr1,chr2,...,chr22`                 | List of chromosomes (e.g., `chr1,chr2,chr3,...`) |
| `ichorcna_readcounter_quality`      | Minimum read count quality to keep in ichorCNA                              | `20`                                  | Any integer                                      |
| `ichorcna_chrs_to_use`              | Chromosomes to use for ichorCNA                                             | `paste0('chr', c(1:22))`              | String format (e.g., `"chr1,chr2,...,chr22"`)    |
| `ichorcna_chrs_to_train`            | Chromosomes to use for training during an ichorCNA run                      | `paste0('chr', c(1:22))`              | String format (e.g., `"chr1,chr2,...,chr22"`)    |
| `ichorcna_chrs_to_normalize`        | Chromosomes to use for normalization during an ichorCNA run                 | `paste0('chr', c(1:22))`              | String format (e.g., `"chr1,chr2,...,chr22"`)    |
| `ichorcna_estimate_normal`          | Whether ichorCNA should estimate normal contamination or not                | `true`                                | `true`, `false`                                  |
| `ichorcna_fraction_reads_male`      | Fraction of data used for copy number calling                               | `0.001`                               | Any value between 0 and 1                        |
| `ichorcna_male_chrX_logR`           | LogR value for male chromosome X                                            | `0.3`                                 | Any numerical value                              |
| `ichorcna_min_map_score`            | Minimum mapping score for reads to be used by ichorCNA                      | `0.75`                                | Any value between 0 and 1                        |
| `ichorcna_max_frac_genome_subclone` | Maximum fraction of genome allowed for subclone                             | `0.5`                                 | Any value between 0 and 1                        |
| `ichorcna_max_frac_cna_subclone`    | Maximum fraction of CNA allowed for subclone                                | `0.7`                                 | Any value between 0 and 1                        |
| `ichorcna_min_segment_bins`         | Minimum number of bins required for segmenting                              | `50`                                  | Any positive integer                             |
| `ichorcna_max_cn`                   | Maximum copy number to be considered by ichorCNA                            | `5`                                   | Any positive integer                             |
| `ichorcna_include_homd`             | Call also homozygous deletions in ichorCNA                                  | `"FALSE"`                             | `"TRUE"`, `"FALSE"`                              |
| `ichorcna_txne`                     | Tumor-normal estimation strength for ichorCNA                               | `0.9999`                              | Any value between 0 and 1                        |
| `ichorcna_alt_frac_threshold`       | Alternative fraction threshold for ichorCNA                                 | `0.05`                                | Any value between 0 and 1                        |
| `ichorcna_trx_strength`             | Strength of transcription effect for ichorCNA                               | `10000`                               | Any integer                                      |
| `ichorcna_plotfiletype`             | Output plot file type for ichorCNA                                          | `pdf`                                 | `pdf`, `png`, `jpeg`, etc.                       |
| `ichorcna_plotylim`                 | Y-axis limits for the generated plots                                       | `c(-2,4)`                             | String format (e.g., `"c(-2, 4)"`)               |
| `ichorcna_estimate_sc`              | Estimate copy number subclonality in ichorCNA                               | `false`                               | `true`, `false`                                  |
| `ichorcna_estimate_ploidy`          | Estimate ploidy in ichorCNA                                                 | `true`                                | `true`, `false`                                  |
| `ichorcna_filter_bam_pon`           | Apply PON filtering on BAM files                                            | N/A                                   | `true`, `false`                                  |
| `ichorcna_normal_states`            | Fraction of normal copy number states used in ichorCNA                      | `0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99` | List of numerical values (e.g., `0.5, 0.6, 0.7`) |
| `ichorcna_gc_wig`                   | Path to a Wiggle file with GC content data for the specified genome         | N/A                                   | File path (must be a `.wig` file)                |
| `ichorcna_map_wig`                  | Path to a Wiggle file with mappability scores for the specified genome      | N/A                                   | File path (must be a `.wig` file)                |
| `ichorcna_reptime_wig`              | Path to a Wiggle file with replication timing data for the specified genome | N/A                                   | File path (must be a `.wig` file)                |
| `ichorcna_centromere_file`          | Path to a file with centromere data for the specified genome                | N/A                                   | File path (must be a `.txt` or similar format)   |

#### _WisecondorX specific options_

| Parameter Name          | Description                                                                      | Default Value | Possible Values                      |
| ----------------------- | -------------------------------------------------------------------------------- | ------------- | ------------------------------------ |
| `wisecondorx_blacklist` | Path to a BED file including loci not to be included in the WisecondorX analysis |               | File path (must be a valid BED file) |
| `wisecondorx_no_rm_dup` | Don't remove duplicates from BAM files                                           |               | `true`, `false`                      |
| `wisecondorx_yfrac`     | Fraction of data used for copy number calling                                    | 0.4           | Any value between 0 and 1            |
| `wisecondorx_zscore`    | Z-score threshold for copy number calling                                        | 5             | Any positive integer                 |
| `wisecondorx_ylim`      | Y axis limits for the generated plots                                            |               | String format (e.g., "0,100")        |

## Downstram analysis options

### Copy Number Signatures

> _NB: CIN Signatures analysis can be performed only with with `--caller ascat_sc` parameter_

| Parameter Name       | Description                                                                                       | Possible Values |
| -------------------- | ------------------------------------------------------------------------------------------------- | --------------- |
| `compute_signatures` | Compute [chromosomal instability signatures](https://www.nature.com/articles/s41586-022-04789-9). | `true`, `false` |

### GISTIC2.0 Analysis

> _NB: Currently, GISTIC analysis can be performed only with genome builds `--genome hg38` or `--genome hg19`_

| Parameter Name            | Description                                                      | Default Value | Possible Values                   |
| ------------------------- | ---------------------------------------------------------------- | ------------- | --------------------------------- |
| `run_gistic`              | Run GISTIC analysis.                                             |               | `true`, `false`                   |
| `gistic_t_amp`            | Default log2ratio threshold to call amplifications.              | 0.1           | Any number (log2 ratio threshold) |
| `gistic_t_del`            | Default log2ratio threshold to call deletions.                   | 0.1           | Any number (log2 ratio threshold) |
| `gistic_remove_x`         | Whether to remove or not chromosome X from the analysis.         |               | `true`, `false`                   |
| `gistic_conf`             | GISTIC confidence level for calling recurrent altered regions.   | 0.99          | Any number (confidence level)     |
| `gistic_qval`             | Maximum q-value to call a region significant.                    | 0.05          | Any number (q-value threshold)    |
| `gistic_broad_analysis`   | Run arm-level (broad) analysis.                                  |               | `true`, `false`                   |
| `gistic_broad_chr_length` | Fraction of altered chromosome to be included in broad analysis. | 0.99          | Any number (fraction)             |
| `gistic_cn_cap`           | GISTIC maximum copy number (higher values will be floored).      | 6             | Any integer (copy number)         |

## nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
