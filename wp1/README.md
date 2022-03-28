# Work Package 1

## Introduction

WP1 provides a pipeline to preprocess sequencing data of different tissues. Originally, the direct sequencing data were to be used for this purpose. However, this pipeline starts after alignment and peak calling of the different reads. Thus, .h5ad files are assumed and the following steps are performed:
- Filtering of cells (by min/max features, features without associated cells, etc.).
- Clustering
- Gene Annotation
- Cell Type Annotation (incl. filtering and ranking)

The most of preprocessing steps has been done with (Epi-)Scanpy and contains Filtering, Clustering and Gene Annotation. Cell Type Annotation was identified with SCSA which brings its own database. Afterwards, the results are written to different files and made available for the other Work Packages.

## Requirements

- Anaconda 4.10.3
- Python 3.9.0
- Relevant packages (also listed in wp1/preprocessing.yml)
    - ipython 8.0.1
    - scanpy 1.8.0
    - episcanpy 0.3.2
    - scikit-learn 1.0.2
    - pandas 1.3.5
    - numpy 1.21.5
    - anndata 0.7.8
- SCSA 1.0
- ~64GB of RAM (recommended)

## Setup

- First set up a new environment with anaconda by the provided `wp1/preprocessing.yml` using the following command
```Bash
conda env create -n "[name_of_your_env]" -f wp1/preprocessing.yml
```
- Next add your new environment as a Jupyter Notebook kernel with
```Bash
python -m ipykernel install --user --name [name_of_your_env] --display-name "[displayed_name]"
```
- At last start a new Jupyter Notebook

## Input Files

- tissue.h5ad files of tissues (A tissue can be divided into smaller samples)
- homo_sapiens.104.main.Chr.gtf

## Output

- .csv and .mtx files about further tissue informations -> WP4
- .bed files (with tissue, cluster & celltype informations) -> WP6

## Example Case

An example case is given in the provided `wp1/complete_preprocessing.ipynb` with comments to most of the executed steps.
