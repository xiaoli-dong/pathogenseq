
<!-- # ![nf-core/pathogen](docs/images/nf-core-pathogen_logo.png) -->

<!-- [![GitHub Actions CI Status](https://github.com/nf-core/pathogen/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/pathogen/actions?query=workflow%3A%22nf-core+CI%22)
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
-->
## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->
**pathogenseq** is a pathogen whole genome sequence (WGS) data analysis pipeline, which inclues sequence quality checking, quality control, taxonomy assignment, assembly, assembly quality assessment, assembled contig annotation, mlst, antimicrobial resistance, virulome, plasmid, and taxonomy prediciton.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. <!-- The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community! -->

<p align="center">
<img src="https://github.com/xiaoli-dong/pathogenseq/assets/52679027/e80d13cc-c9a8-46a0-add0-8dd5e23729b7">
</p>
## Pipeline summary

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->
By default, the pipeline supports both short and long reads:

- Sequence quality check and quality control
  - Short reads
    - Short Illumina reads quality checks ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
    - Short read quality control ([BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) | [fastp](https://github.com/OpenGene/fastp))
    - Short read statistics ([seqkit stats](https://bioinf.shenwei.me/seqkit/usage/#stats))
    - Taxonomic assignment and contamination check ([`Kraken2`](https://ccb.jhu.edu/software/kraken2/))
  - Long reads
    - Nanopore long read quality checks ([NanoPlot](https://github.com/wdecoster/NanoPlot))
    - Nanopore long read adapter trimming ([Porechop](https://github.com/rrwick/Porechop))
    - Nanopore long read quality and length filter ([chopper](https://github.com/wdecoster/chopper))
    - Nanopore long read statistics ([seqkit stats](https://bioinf.shenwei.me/seqkit/usage/#stats))
- Assembly
  - Short read assembly with user choice of the assemblers ([Spades](https://github.com/ablab/spades) | [Skesa](https://github.com/ncbi/SKESA) | [Unicycler](https://github.com/rrwick/Unicycler) | [megahit](https://github.com/voutcn/megahit))
  - Long read assembly is following the steps below:
    - Nanopore long read de novo assembly ([`Flye`](https://github.com/fenderglass/Flye))
    - Circular Flye contigs are rotated to start in the center of the contig ([in-house perl script](https://github.com/xiaoli-dong/pathogenseq/blob/main/bin/reset_start_position_for_circular_genome.pl))
    - Long read polishing and consensus generating ([`Medaka`](https://github.com/nanoporetech/medaka))
    - Short-read polishing while short reads are available:
      - [Polypolish](https://github.com/rrwick/Polypolish)
      - [POLCA](https://github.com/alekseyzimin/masurca)
- Assembly quality check
  - Rapid assessment of genome assembly completeness and contamination using machine learning approach ([CheckM2](https://github.com/chklovski/CheckM2))
  - Rapid taxonomic identification of microbial pathogens from assemblies and also the assement of the sample relatedness ([gambit](https://github.com/jlumpe/gambit))
- Genome annotation
  - Gene prediction and annotation ([Bakta](https://github.com/oschwengers/bakta))
  - Identify acquired antimicrobial resistance genes in the assembled contigs ([AMRFinderPlus](https://github.com/ncbi/amr))
  - Scan contig files against traditional PubMLST typing schemes ([mlst](https://github.com/tseemann/mlst))
  - Typing and reconstruction of plasmid sequences from assembled contigs ([MOB-suite](https://github.com/phac-nml/mob-suite))
  - Virulome detection ([abricate](https://github.com/tseemann/abricate) with [VFDB](http://www.mgc.ac.cn/VFs/main.htm))
- Summarize and generate the analysis report, software version control reports

## Pipeline reference databases
* [Bakta database](https://github.com/oschwengers/bakta#database)
* [Kraken2 database](https://benlangmead.github.io/aws-indexes/k2)
* [AMRFinderPlus database](https://github.com/ncbi/amr/wiki/AMRFinderPlus-database)
* [CheckM2 database](https://github.com/chklovski/CheckM2)
* [gambit database](https://github.com/jlumpe/gambit)
  
## Quick Start

The workflow uses nextflow to manage compute and software resources, as such nextflow will need to be installed before attempting to run the workflow.
The workflow can currently be run using either singularity or conda to provide isolation of the required software. Both methods are automated out-of-the-box provided either docker or singularity is installed.
It is not required to clone or download the git repository in order to run the workflow. 

  > **Note**
  > If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
  > to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
  > with `-profile test` before running the workflow on actual data.
 
  > Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_



### Check workflow options
You can clone or download the pathogenseq from github to local computer or you can directly run the pipeline from github. To check the pipeline command line options:

```{r df-drop-ok, class.source="bg-success"}
# running directly from github without downloading or cloning
nextflow run xiaoli-dong/pathogenseq -r revision_number(e.g:8657a20) --help

# download the pipeine and run the program from the local computer
nextflow run your_path_to/pathogenseq/main.nf --help
N E X T F L O W  ~  version 23.04.1
Launching `main.nf` [happy_linnaeus] DSL2 - revision: 4ef093f544
WARN: Access to undefined parameter `show_hidden_params` -- Initialise it to a default value eg. `params.show_hidden_params = some_value`


------------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
  nf-core/pathogenseq v1.0.1
------------------------------------------------------
Typical pipeline command:

  nextflow run nf-core/pathogenseq --input samplesheet.csv --outdir results --platform illumina -profile singularity

Input/output options
  --input                        [string]  Path to comma-separated samplesheet file containing information about the samples in the experiment.
  --outdir                       [string]  The output directory where the results will be saved. You have to use absolute paths to storage on Cloud 
                                           infrastructure. 
  --platform                     [string]  Specifies the platform used to generate the sequences - available options are 'illumina|nanopore'. [default: 
                                           illumina] 
  --igenomes_ignore              [boolean] Whether ignore igenome configuration loading.
  --email                        [string]  Email address for completion summary.
  --multiqc_title                [string]  MultiQC report title. Printed as page header, used for filename if not otherwise specified.

illumina_options
  --illumina_reads_qc_tool       [string]  Specifies the reads triming and qc tool to use - available options are 'fastp|bbduk'. [default: fastp]
  --illumina_reads_assembler     [string]  Specifies the illumina reads assembly tool to use - available options are 'megahit|spades|skesa|unicycler'. 
                                           [default: unicycler] 
  --skip_illumina_reads_qc       [boolean] Skip illumina read quality control step. [default: false]
  --skip_illumina_reads_assembly [boolean] Skip illumina read assembly step. [default: false]
  --skip_illumina_kraken2        [boolean] Skip kraken2 with illumina reads step. [default: false]

nanopore_options
  --skip_nanopore_kraken2        [boolean] Skip kraken2 with nanopore reads step. [default: true]
  --nanopore_reads_assembler     [string]  Specifies the nanopore reads assembly tool to use - available options are 'flye+medaka'. [default: 
                                           flye+medaka] 
  --skip_polypolish              [boolean] Skip assembly polish step with illumina reads with polypolish tool. [default: false]
  --skip_polca                   [boolean] Skip assembly polish step with illumina reads with polca tool. [default: false]
  --skip_nanopore_reads_qc       [boolean] Skip nanopore read quality control step. [default: false]
  --skip_nanopore_reads_assembly [boolean] Skip nanopore read assembly step. [default: false]
  --skip_illumina_reads_polish   [boolean] Skip assembly polish step with illumina reads steps. [default: false]

annotation_options
  --bakta_db                     [string]  Path to bakta database. [default: /nfs/APL_Genomics/db/prod/bakta/db]
  --checkm2_db                   [string]  Path to checkm2 database. [default: /nfs/APL_Genomics/db/prod/CheckM2_database/uniref100.KO.1.dmnd]
  --amrfinderplus_db             [string]  null
  --kraken2_db                   [string]  Specify path to kraken2 database [default: /nfs/APL_Genomics/db/prod/kraken2/k2_standard_08gb_20220926]
  --gambit_db                    [string]  Path to gambit database. [default: /nfs/APL_Genomics/db/prod/gambit]
  --skip_checkm2                 [boolean] Skip checkm2 step. [default: false]
  --skip_gambit                  [boolean] Skip gambit step. [default: false]
  --skip_bakta                   [boolean] Skip bakta step. [default: false]
  --skip_mlst                    [boolean] Skip mlst step. [default: false]
  --skip_mobsuite                [boolean] Skip skip_mobsuite step. [default: false]
  --skip_virulome                [boolean] Skip virulome step. [default: false]
  --skip_multiqc                 [boolean] null [default: false]
  --skip_amr                     [boolean] Skip amr step. [default: false]
  --skip_depth_and_coverage      [boolean] Skip assembly depth calculation step. [default: false]

!! Hiding 22 params, use --show_hidden_params to show them !!
------------------------------------------------------
If you use nf-core/pathogenseq for your analysis please cite:

* The nf-core framework
  https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
  https://github.com/nf-core/pathogenseq/blob/master/CITATIONS.md
------------------------------------------------------
```


### Prepare required samplesheet input
The pathogenseq pipeline requires user to provide a csv format samplesheet, which contains the sequenence information for each sample, as input. See below for what the samplesheet looks like:

```samplesheet.csv```

```
sample,fastq_1,fastq_2,long_fastq,basecaller_mode
sample1,shortreads_1.fastq.gz,shortreads_2.fastq.gz,longreads.fastq.gz,r1041_e82_400bps_hac_v4.2.0
sample2,shortreads.fastq,NA,longreads.fastq.gz,r1041_e82_400bps_sup_v4.2.0
sample3,NA,NA,longreads.fastq.gz,NA
sample4,shortreads_1.fastq.gz,shortreads_2.fastq.gz,NA
```
The csv format samplesheet has five required columns:
* The first row of the csv file is the header describing the columns
* Each row represents a unique sample to be processed, the first colum is the unique sample id
* When the information for a particular column is missing, please fill the column with "NA"
* The "fastq_1" and "fastq_2" columns are reserved for supplying the short sequence files
* "basecaller_mode" is for user to provide the Nanopore basecalling model, for example: r1041_e82_400bps_hac_v4.2.0. The availble models for medaka are listed at bewlow:

```
xdong@M691822:~$ medaka tools list_models
Available: r103_fast_g507, r103_fast_snp_g507, r103_fast_variant_g507, r103_hac_g507, r103_hac_snp_g507, r103_hac_variant_g507, r103_min_high_g345, r103_min_high_g360, r103_prom_high_g360, r103_prom_snp_g3210, r103_prom_variant_g3210, r103_sup_g507, r103_sup_snp_g507, r103_sup_variant_g507, r1041_e82_260bps_fast_g632, r1041_e82_260bps_fast_variant_g632, r1041_e82_260bps_hac_g632, r1041_e82_260bps_hac_v4.0.0, r1041_e82_260bps_hac_v4.1.0, r1041_e82_260bps_hac_variant_g632, r1041_e82_260bps_hac_variant_v4.1.0, r1041_e82_260bps_sup_g632, r1041_e82_260bps_sup_v4.0.0, r1041_e82_260bps_sup_v4.1.0, r1041_e82_260bps_sup_variant_g632, r1041_e82_260bps_sup_variant_v4.1.0, r1041_e82_400bps_fast_g615, r1041_e82_400bps_fast_g632, r1041_e82_400bps_fast_variant_g615, r1041_e82_400bps_fast_variant_g632, r1041_e82_400bps_hac_g615, r1041_e82_400bps_hac_g632, r1041_e82_400bps_hac_v4.0.0, r1041_e82_400bps_hac_v4.1.0, r1041_e82_400bps_hac_v4.2.0, r1041_e82_400bps_hac_variant_g615, r1041_e82_400bps_hac_variant_g632, r1041_e82_400bps_hac_variant_v4.1.0, r1041_e82_400bps_hac_variant_v4.2.0, r1041_e82_400bps_sup_g615, r1041_e82_400bps_sup_v4.0.0, r1041_e82_400bps_sup_v4.1.0, r1041_e82_400bps_sup_v4.2.0, r1041_e82_400bps_sup_variant_g615, r1041_e82_400bps_sup_variant_v4.1.0, r1041_e82_400bps_sup_variant_v4.2.0, r104_e81_fast_g5015, r104_e81_fast_variant_g5015, r104_e81_hac_g5015, r104_e81_hac_variant_g5015, r104_e81_sup_g5015, r104_e81_sup_g610, r104_e81_sup_variant_g610, r10_min_high_g303, r10_min_high_g340, r941_e81_fast_g514, r941_e81_fast_variant_g514, r941_e81_hac_g514, r941_e81_hac_variant_g514, r941_e81_sup_g514, r941_e81_sup_variant_g514, r941_min_fast_g303, r941_min_fast_g507, r941_min_fast_snp_g507, r941_min_fast_variant_g507, r941_min_hac_g507, r941_min_hac_snp_g507, r941_min_hac_variant_g507, r941_min_high_g303, r941_min_high_g330, r941_min_high_g340_rle, r941_min_high_g344, r941_min_high_g351, r941_min_high_g360, r941_min_sup_g507, r941_min_sup_snp_g507, r941_min_sup_variant_g507, r941_prom_fast_g303, r941_prom_fast_g507, r941_prom_fast_snp_g507, r941_prom_fast_variant_g507, r941_prom_hac_g507, r941_prom_hac_snp_g507, r941_prom_hac_variant_g507, r941_prom_high_g303, r941_prom_high_g330, r941_prom_high_g344, r941_prom_high_g360, r941_prom_high_g4011, r941_prom_snp_g303, r941_prom_snp_g322, r941_prom_snp_g360, r941_prom_sup_g507, r941_prom_sup_snp_g507, r941_prom_sup_variant_g507, r941_prom_variant_g303, r941_prom_variant_g322, r941_prom_variant_g360, r941_sup_plant_g610, r941_sup_plant_variant_g610
Default consensus:  r1041_e82_400bps_sup_v4.2.0
Default variant:  r1041_e82_400bps_sup_variant_v4.2.0
```

### Run the pipeline:
> * If you are using `singularity` then the pipeline will auto-detect this and attempt to download the Singularity images directly as opposed to performing a conversion from Docker images. If you are persistently observing issues downloading Singularity images directly due to timeout or network issues then please use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. Alternatively, it is highly recommended to use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to pre-download all of the required containers before running the pipeline and to set the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options to be able to store and re-use the images from a central location for future pipeline runs.
> * If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

```
# Runing the pipeline from remote github
nextflow run xiaoli-dong/pathogenseq \
  -r seven_character_github_revision_number (e.g: 8657a20)
  --input samplesheet.csv \
  -profile <docker|singularity|podman|shifter|charliecloud|conda/institute> \
  --outdir results \
  --platform <illumina|nanopore>

# run the pipeline with a local clone
nextflow run your_path_to/pathogenseq/main.nf \
  --input samplesheet.csv \
  -profile <docker|singularity|podman|shifter|charliecloud|conda/institute> \
  --outdir results \
  --platform <illumina|nanopore>

# an example command to launch the pipeline from local computer and run it with ```singularity``` configraton profile.
nextflow run your_path_to/pathogenseq/main.nf \
  -profile singularity \
  --input samplesheet.csv \
  --outdir results_singularity \
  --platform <illumina|nanopore>

# an example commnad to launch the pipeline from a local clone and run it with ```conda``` configraton profile.
nextflow run your_path_to/pathogenseq/main.nf \
  -profile conda \
  --input samplesheet.csv \
  --outdir results_conda \
  --platform <illumina|nanopore>
```
>* Notes: Please provide pipeline parameters via the CLI or Nextflow -params-file option. Custom config files including those provided by the -c Nextflow option can be used to provide any configuration except for parameters; see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).


## Pipeline output
<table>
  <tr valign="top">
    <td>
      <p>Fig 1. Pathogenseq output top level layout</p>
      <img src="https://github.com/xiaoli-dong/pathogenseq/assets/52679027/5e6d5e43-fda1-43bc-b5c4-d844e8f1a7cf">
    </td>
    <td>
       <p>Fig 2. Pathogenseq pipeline_info directory layout</p>
     <img src="https://github.com/xiaoli-dong/pathogenseq/assets/52679027/1ccbb832-6b97-4ab4-bde1-d5ecdebd3fbe">
    </td>
    <td>
       <p>Fig 3. Pathogenseq report directory layout</p>
     <img src="https://github.com/xiaoli-dong/pathogenseq/assets/52679027/59ecc195-2b3c-4b78-951b-38e7595c85fd">
    </td>
   
  </tr>
</table>
From the above screenshots, you can see:

* Results of each sample included in the analysis go to its own directry. 
* pipeline_info directory contains software version control information, nextflow workflow report, resource usage report, task report. 
* report directory contains the summary files of each analysis task

<table >
  <caption>Example pathogenseq data analysis outputs for a particalur sample</caption>
  <tr valign="top">
    <td>
     <img src="https://github.com/xiaoli-dong/pathogenseq/assets/52679027/ace1cbc3-ccd3-4cf3-a2ef-dc36c9651bf2" width="120%">
    </td>
    <td>
     <img src="https://github.com/xiaoli-dong/pathogenseq/assets/52679027/2422bf48-3089-4e68-8472-5573c2329929" width="120%">
    </td>
  </tr>
</table>
 

## Credits
pathogenseq was originally written by Xiaoli Dong. Extensive support was provided from other co-authors on the scientific or technical input required for the pipeline:
* Dr. Matthew Croxen
* Dr. Tarah Lynch

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

<!-- 
## Citations
-->
<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/pathogen for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->
<!-- An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
-->
## Notes
with conda, fastqc 0.11.9--0 is not working. I need to change the version to 0.12.1. 

mobsuite is not working with 3.0.3 in conda, I need to change the mobsuite version to  3.1.4
## Reference
<!-- >>https://github.com/ablab/graphamr -->
