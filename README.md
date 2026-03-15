[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![](https://img.shields.io/static/v1?label=CLI&message=Snaketool&color=blueviolet)](https://github.com/beardymcjohnface/Snaketool)
![GitHub last commit (branch)](https://img.shields.io/github/last-commit/npbhavya/MitoBee?branch=main)

<img src="test-files/mitobeelogo.png" width="200" alt="MitoBee logo">

# MitoBee

## Snakemake workflow to get mitogenomes from metagenomic data
Still under development! Stable release out as a version, but only if there is a closely related mitogenome available.

Documentation: [Wiki](https://github.com/npbhavya/Mitobee/wiki/Home)

### Install 

**Source install**
Run the below commands:

    git clone https://github.com/npbhavya/MitoBee.git
    cd MitBee
    mamba create -y -n mitobee python=3.13
    conda activate mitobee
    pip install -e . 

**Once I have a stable version release, I will upload them to conda and pip as well**

### Running the code
Note: This code works only on paired end metagenomes for now. 

This workflow is made modular; 

1. `mitobee run` 
If there is a representative closely related genome mitogenome, provide that as the host seq and get started

```
    #Running mitobee with test files available in the repo
    mitobee run --input test-files/metagenomes --extn fastq.gz \
         --pattern_r1 _R1 --pattern_r2 _R2 \
         --host_seq test-files/am-dh4.fasta \
         --output output
```

2. `mitobee tree`
Once the mitochondrial genomes are built from each metagenome sample, run this module to build a tree with these mitogenomes and other references

```
    #After the mitogenomes are made from the mitobee run results. Add other references to build a tree
    #Once again example with test files
    mitobee tree --input test-files/mitogenomes --extn fasta --output output -k all

```

3.  `mitobee search`
If there is a no closely related mitogenome available, then this step can be run first to search against a set of mitogenomes or mito genes 
This module will provide an overview of which reference to use

```
    #If a closely related mitochondrial genome is not available, but a gene is, like cox  or rRNA genes
    #Download the reference genes you would like to use of the closely related genomes

    #to search against mitogenomes refernece set
    mitobee search --input test-files/mitogenomes --extn fasta --ref_seq  test-files/ref --output output -k all --mode mitogenome

    #to search against mitogenomes refernece gene set 
    mitobee search --input test-files/mitogenomes --extn fasta --ref_seq  test-files/ref --output output -k all --mode genes 

```

### Input files
Input files:
- Input directory with metagenomes
- Reference directory
  - If running `run`  or `tree` module, provide a **(one)** reference genome. 
  - If running `gene` module, provide a reference gene set

### Output files
Output files: Provide the output folder, contains subdirectories
- PROCESSING: Folder containing intermediate files
- REPORTS: Final results including the mitogenome fasta files from (hopefully) each metagenome sample \
      Also inlcudes the QC reports, to include stats on how many reads were processed, and not
