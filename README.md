# Morphmap

## Table of contents

- [Metadata](#metadata)
- [Images](#images)
    * [How to download images?](#how-to-download-images)
- [Extracting single-cell features using CellProfiler](#extracting-single-cell-features-using-cellprofiler)
    * [Creation of single-cell profiles](#creation-of-single-cell-profiles)
        * [How to download single-cell profiles?](#how-to-download-single-cell-profiles) 
    * [Creation of well-level profiles](#creation-of-well-level-profiles)
        * [How to download well-level profiles?](#how-to-download-well-level-profiles)
- [Processing the profiles](#processing-the-profiles)
- [How to run the analyses?](#how-to-run-the-analyses)
    * [Cloning this repo](#cloning-this-repo)
    * [Downloading the profiles](#downloading-the-profiles)
    * [Running the analyses notebooks](#running-the-analyses-notebooks)
        * [Installing the conda environment](#installing-the-conda-environment)


<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>

**Link to the biorxiv manuscript: https://www.biorxiv.org/content/10.1101/2024.12.02.624527 **

This is a dataset of images and profiles generated as a part of the [JUMP Cell Painting (JUMP-CP) project](https://jump-cellpainting.broadinstitute.org/). Genes were either over-expressed (ORF) or knocked out (CRISPR) and the cells were assayed using an imaging assay called [Cell Painting](https://jump-cellpainting.broadinstitute.org/cell-painting). From the images, features were extracted using the [CellProfiler](https://cellprofiler.org/) software. The features were then processed and the resulting profiles were the analyzed using notebooks in this repository. 

**Many resources and additional views / formats of the data are available at the [JUMP Cell Painting Hub](https://broad.io/jump).**

In the following sections, instructions are provided for downloading the various components of this dataset, processing the dataset and analyzing the profiles.

## Metadata
Metadata information, such as, which plate from which batch contains a particular gene, is available in the [datasets repo](https://github.com/jump-cellpainting/datasets/tree/main/metadata).

## Images
Cell images are available for download from the [cellpainting gallery](https://cellpainting-gallery.s3.amazonaws.com/index.html) public AWS S3 bucket. 

### How to download images?
There are two `sources` of data. ORF images are from `source_4` and CRISPR images are from `source_13`.

```bash
source=<SOURCE NAME>
aws s3 sync \
  --no-sign-request \
  s3://cellpainting-gallery/cpg0016-jump/${source}/images/ . 
```

## Extracting single-cell features using CellProfiler
Features were extracted using the CellProfiler pipeline in https://github.com/broadinstitute/imaging-platform-pipelines/tree/master/JUMP_production#production-pipelines. 

### Creation of single-cell profiles
Instructions for creating the single-cell profiles from images are provided in the [Image-based profiling handbook](https://cytomining.github.io/profiling-handbook/01-overview.html).

#### How to download single-cell profiles? 
Single-cell profiles can be downloaded from the [cellpainting gallery](https://cellpainting-gallery.s3.amazonaws.com/index.html) public AWS S3 bucket.

```bash
source=<SOURCE NAME>
batch=<BATCH NAME>
plate=<PLATE NAME>
aws s3 sync \
  --no-sign-request \
  s3://cellpainting-gallery/cpg0016-jump/${source}/workspace/backend/${batch}/${plate}/ --exclude "*" --include "*.sqlite" .
```

### Creation of well-level profiles
Well-level profiles are also created using the instructions provided in the [Image-based profiling handbook](https://cytomining.github.io/profiling-handbook/01-overview.html).

#### How to download well-level profiles?
Well-level profiles can also be downloaded from the [cellpainting gallery](https://cellpainting-gallery.s3.amazonaws.com/index.html) public AWS S3 bucket.

```bash
source=<SOURCE NAME>
batch=<BATCH NAME>
plate=<PLATE NAME>
aws s3 sync \
  --no-sign-request \
  s3://cellpainting-gallery/cpg0016-jump/${source}/workspace/profiles/${batch}/${plate}/ .
```

## Processing the profiles
Various steps were performed to remove technical noise from the profiles. These steps are as follows:

- Well position correction
- Cell count regression
- Normalization
- Outlier removal
- Feature selection
- Sphering
- Harmony correction

These steps can be performed using the [jump-profiling-recipe](https://github.com/broadinstitute/jump-profiling-recipe/tree/v0.1.0) and the appropriate config file (`orf.json` and `crispr.json` from the `input` folder).

The processed profiles are stored in the [cellpainting-gallery bucket](https://cellpainting-gallery.s3.amazonaws.com/index.html#cpg0016-jump-assembled/source_all/workspace/profiles/jump-profiling-recipe_2024_a917fa7/). 

## How to run the analyses?
### Cloning this repo
To download/clone this repository, run the following commands

```bash
git clone https://github.com/jump-cellpainting/2024_Chandrasekaran_Morphmap.git 
cd 2024_Chandrasekaran_Morphmap
git submodule update --init --recursive
```

### Downloading the profiles
To download the profiles and other files required to run the analyses in this repository, run the following commands

```bash
cd profiles
./download-profiles.sh
```

### Running the analyses notebooks
Notebooks in each folder are run to replicate the analyses. To reproduce the analyses, install the conda environment within each folder and run the notebooks.

#### Installing the conda environment
Download and install mamba from [miniforge](https://github.com/conda-forge/miniforge) using the appropriate installer for your operating system.

Once installed, run the following command to create the conda environment in each folder using the following commands

```bash
mamba env create -f environment.yml
mamba activate <conda environment name>
```
