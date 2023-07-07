<!-- # ![nf-core/pathogen](docs/images/nf-core-pathogen_logo.png) -->

[![GitHub Actions CI Status](https://github.com/nf-core/pathogen/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/pathogen/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/pathogen/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/pathogen/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/pathogen/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23pathogen-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/pathogen)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->
**nf-core/pathogenseq** is a bioinformatics best-practice analysis pipeline for Pathogen sequence quality checking, quality control, assembly, annotation, predicting mlst, antimicrobial resistance, virulome, and plasmid.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!


## Pipeline summary

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->
The pipeline currently performs the following:

* Sequence quality check and quality control
  * Short reads
    * Short Illumina reads quality checks ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
    * Short read quality control ([BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) | [fastp](https://github.com/OpenGene/fastp))
    * Short read statistics ([seqkit stats](https://bioinf.shenwei.me/seqkit/usage/#stats))
    * Taxonomic assignment and contamination check ([`Kraken2`](https://ccb.jhu.edu/software/kraken2/))
  * Long reads
    * Nanopore long read quality checks ([NanoPlot](https://github.com/wdecoster/NanoPlot))
    * Nanopore long read adapter trimming ([Porechop](https://github.com/rrwick/Porechop))
    * Nanopore long read statistics ([seqkit stats](https://bioinf.shenwei.me/seqkit/usage/#stats))
* Assembly
  * Short read assembly with user choice of the assemblers ([Spades](https://github.com/ablab/spades) | [Skesa](https://github.com/ncbi/SKESA) | [Unicycler](https://github.com/rrwick/Unicycler) | [megahit](https://github.com/voutcn/megahit) | [shovill](https://github.com/tseemann/shovill) | [shovill_skesa](https://github.com/tseemann/shovill) | [shovill_megahit](https://github.com/tseemann/shovill) | [shovill_skesa](https://github.com/tseemann/shovill)
  * Long read assembly is following the steps below:
    * Nanopore long read de novo assembly ([`Flye`](https://github.com/fenderglass/Flye))
    * Circular Flye contigs are rotated to start in the center of the contig ([in-house perl script](https://github.com/xiaoli-dong/pathogenseq/blob/main/bin/reset_start_position_for_circular_genome.pl))
    * Long read polishing and consensus generating ([`Medaka`](https://github.com/nanoporetech/medaka))
  * Short-read polishing if short read available: 
    * [Polypolish](https://github.com/rrwick/Polypolish)
    * [POLCA](https://github.com/alekseyzimin/masurca)
      
* Contig annotation [Bakta](https://github.com/oschwengers/bakta)
* Identify acquired antimicrobial resistance genes in the assembled contigs ([AMRFinderPlus](https://github.com/ncbi/amr))
* Scan contig files against traditional PubMLST typing schemes [mlst](https://github.com/tseemann/mlst)
* Typing and reconstruction of plasmid sequences from assembled contigs [MOB-suite](https://github.com/phac-nml/mob-suite)
* Virulome detection ([abricate](https://github.com/tseemann/abricate) with [VFDB](http://www.mgc.ac.cn/VFs/main.htm))
* Summarize and generate the analysis report, software version control report

## Usage

 If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with -profile test before running the workflow on actual data.


## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```console
    nextflow run xiaoli-dong/pathogen -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

    > * Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
    > * If you are using `singularity` then the pipeline will auto-detect this and attempt to download the Singularity images directly as opposed to performing a conversion from Docker images. If you are persistently observing issues downloading Singularity images directly due to timeout or network issues then please use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. Alternatively, it is highly recommended to use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to pre-download all of the required containers before running the pipeline and to set the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options to be able to store and re-use the images from a central location for future pipeline runs.
    > * If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

    <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

    ```console
    nextflow run xiaoli-dong/pathogen -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> -c conf/example.config -profile singularity
    ```

## Documentation

The nf-core/pathogen pipeline comes with documentation about the pipeline [usage](https://nf-co.re/pathogen/usage), [parameters](https://nf-co.re/pathogen/parameters) and [output](https://nf-co.re/pathogen/output).

## Credits

nf-core/pathogen was originally written by Xiaoli Dong.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#pathogen` channel](https://nfcore.slack.com/channels/pathogen) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/pathogen for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->
An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).


## bugs:

> POLCA is not working with sigularity, bwa is not the path
> nextpolish has no singularity image, need local installation
>need Tiny in the system perl

## reference
<!-- >>https://github.com/ablab/graphamr -->
