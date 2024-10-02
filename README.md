# SAMURAI: Shallow Analysis of Copy nuMber alterations Using a Reproducible And Integrated bioinformatics pipeline <img src="https://github.com/DIncalciLab/samurai/assets/71792548/7e842710-6fa6-4f0e-997d-33f95693158a" width="150" height="150" align="right">

<quote>
戦の勝負は、将と士卒の気によるのみなり。

(The outcome of the battle depends on the spirit of the commander and the men.)

Minamoto no Yoshitsune (1159-1189)
</quote>

## Introduction

**SAMURAI** is a bioinformatics best-practice analysis pipeline for the analysis of shallow whole genome sequencing (sWGS) data for the identification of copy number alterations (CNAs). It supports a number of workflows depending on the nature of the samples (coming from tissues or other biological fluids like plasma). While it was developed with cancer studies in mind, it is applicable to any field where DNA alterations need to be studied.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Pipeline summary

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. (Optional) Trim reaads based on quality scores and extract unique molecular identifiers (UMIs) if applicable ([`fastp`]())
3. Align reads to the reference genome ([`BWA`]())
4. Run quality control checks on the aligned reads ([`Picard`]())
5. Perform copy number alteration identification for tissue samples([`QDNAseq`](), [`ASCAT.sc`]())
6. (Optional) Perform size selection on samples from liquid biopsies([`sambamba`]())
7. Perform copy number alteration identification for liquid biopsy samples([`ichorCNA`](), [`WisecondorX`]())
8. (Optional) Extract copy number instability signatures ([`CINSignatureQuantification`]())
9. (Optional) Identify recurrent altered regions in the sample population ([`GISTIC`]())
10. Present QC for each sep of the pipeline ([`MultiQC`](http://multiqc.info/))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=23.10.1`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```bash
   nextflow run dincalcilab/samurai -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   ```bash
   nextflow run dincalcilab/samurai --input samplesheet.csv --outdir <OUTDIR> --genome GRCh37 -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> <pipeline options>
   ```

## Credits

dincalcilab/samurai was originally written by Sara Potente and Luca Beltrame.

We thank the following people for their extensive assistance in the development of this pipeline:

- Laura Mannarino
- Riccardo Zadro

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

If you use dincalcilab/samurai for your analysis, please cite it using the following doi (currently a preprint): [10.1101/2024.09.30.615766 ](https://doi.org/10.1101/2024.09.30.615766)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
